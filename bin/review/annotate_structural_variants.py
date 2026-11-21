#!/usr/bin/env python3
"""Gene annotation of native SV evidence, with explicit-mate deduplication only."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import pysam
import pysam.bcftools
from breakpoint_annotation import gene_index, overlap


def endpoint(rec):
    alt=rec.alts[0] if rec.alts else ''
    match=re.search(r'([\[\]])([^\[\]]+):(\d+)\1',alt)
    if match:
        return dict(chrom=match[2],pos=int(match[3]),local_side='left' if alt[0] in '[]' else 'right',remote_side='right' if match[1]==']' else 'left',encoding='BND')
    chrom=rec.info['CHR2'] if 'CHR2' in rec.info else None
    pos=rec.info['POS2'] if 'POS2' in rec.info else None
    if isinstance(chrom,tuple):chrom=chrom[0]
    if isinstance(pos,tuple):pos=pos[0]
    strands = rec.info.get('STRANDS')
    if isinstance(strands,tuple):strands=strands[0] if len(strands)==1 else None
    strands=str(strands).split(':')[0] if strands is not None else ''
    local_side='right' if strands[:1]=='+' else 'left' if strands[:1]=='-' else None
    remote_side='left' if strands[1:2]=='+' else 'right' if strands[1:2]=='-' else None
    if chrom is not None and pos is not None:return dict(chrom=str(chrom),pos=int(pos),local_side=local_side,remote_side=remote_side,encoding='CHR2/POS2')
    if rec.stop>rec.pos:return dict(chrom=rec.contig,pos=rec.stop,local_side=local_side,remote_side=remote_side,encoding='interval_END')
    return dict(chrom=None,pos=None,local_side=None,remote_side=None,encoding='unresolved')


def mate_ids(rec):
    for field in ('MATEID','MATE_ID'):
        value=rec.info[field] if field in rec.info else None
        if value:return list(value) if isinstance(value,tuple) else str(value).split(',')
    return []


def annotate(sample,caller,vcf,gff,fai,output,gene_database=None):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    if gene_database:
        import shutil
        shutil.copyfile(gene_database[0],out/'.sv.sqlite')
    db=sqlite3.connect(out/'.sv.sqlite')
    audit=gene_database[1] if gene_database else gene_index(gff,fai,db)
    db.executescript('CREATE TABLE records (ordinal INTEGER PRIMARY KEY,id TEXT,data TEXT,consumed INTEGER DEFAULT 0); CREATE INDEX native_id ON records(id);')
    counts=Counter();review=[];prefix=f'{sample}.{caller}.structural'
    annotated=out/(prefix+'.annotated.vcf.gz')
    with open(vcf,'rb') as native_stream,pysam.VariantFile(native_stream) as source:
        names=list(source.header.samples)
        aliases={sample}
        if caller=='sniffles':aliases.add('SAMPLE')
        if caller in ('severus','stellerator'):aliases.add(sample+'_snp_phased_haplotagged')
        if len(names)!=1 or names[0] not in aliases:raise ValueError('SV sample mismatch: '+str(names))
        native_sample=names[0]
        header=source.header.copy()
        # pysam materializes END when translating symbolic alleles (including
        # Stellerator <BND>). Declare it before the writer serializes its header.
        if 'END' not in header.info:
            header.info.add('END',1,'Integer','End position of the variant described in this record')
        header.info.add('ONCOSEQ_REVIEW_GENES1','.', 'String','Genes overlapping native first breakpoint; no pathogenicity inference')
        header.info.add('ONCOSEQ_REVIEW_GENES2','.', 'String','Genes overlapping native partner breakpoint; no fusion inference')
        with pysam.VariantFile(annotated,'wz',header=header) as writer:
            for ordinal,rec in enumerate(source):
                partner=endpoint(rec);genes1=overlap(db,rec.contig,rec.pos);genes2=overlap(db,partner['chrom'],partner['pos'])
                names1=sorted({g['name'] for g in genes1});names2=sorted({g['name'] for g in genes2})
                native=str(rec).strip();rec.translate(writer.header)
                if names1:rec.info['ONCOSEQ_REVIEW_GENES1']=tuple(names1)
                if names2:rec.info['ONCOSEQ_REVIEW_GENES2']=tuple(names2)
                writer.write(rec)
                interval=[]
                if partner['encoding']=='interval_END':interval=[dict(id=i,name=n) for i,n in db.execute('SELECT id,name FROM genes WHERE chrom=? AND start<=? AND end>=? ORDER BY name,id',(rec.contig,partner['pos'],rec.pos))]
                data=dict(schema_version=1,sample=sample,vcf_sample=native_sample,build='hg38',record_type='structural_junction',id=f'{caller}:{ordinal}',caller=caller,native_id=rec.id,chrom=rec.contig,pos=rec.pos,partner=partner,coordinate_convention='Native VCF 1-based; interval END inclusive',genes1=genes1,genes2=genes2,interval_genes=interval,gene_pair=sorted([names1[0],names2[0]]) if len(names1)==len(names2)==1 else None,mate_ids=mate_ids(rec),filter=list(rec.filter) or ['.'],evidence=dict(info=dict(rec.info),format={k:dict(rec.samples[k]) for k in rec.samples},qual=rec.qual),source_records=[native],clinical_classification=None)
                db.execute('INSERT INTO records(ordinal,id,data) VALUES (?,?,?)',(ordinal,rec.id,json.dumps(data)))
    db.commit()
    # Annotation preserves source record order; sort once before indexing.
    sorted_vcf=out/'.sorted.vcf.gz'
    pysam.bcftools.sort('-Oz','-o',str(sorted_vcf),str(annotated),catch_stdout=False);sorted_vcf.replace(annotated)
    pysam.tabix_index(str(annotated),preset='vcf',force=True)
    with gzip.open(out/(prefix+'.evidence.tsv.gz'),'wt') as stream:
        writer=csv.writer(stream,delimiter='\t');writer.writerow(['id','caller','gene_pair','mate_status','event_json'])
        for ordinal,encoded in db.execute('SELECT ordinal,data FROM records ORDER BY ordinal'):
            if db.execute('SELECT consumed FROM records WHERE ordinal=?',(ordinal,)).fetchone()[0]:continue
            row=json.loads(encoded);row['mate_status']='unpaired'
            if row['mate_ids']:
                candidates=[]
                for mid in row['mate_ids']:candidates+=db.execute('SELECT ordinal,data FROM records WHERE id=? AND consumed=0 AND ordinal!=?',(mid,ordinal)).fetchall()
                if len(candidates)==1:
                    other_ord,other_text=candidates[0];other=json.loads(other_text)
                    reciprocal=(row['native_id'] in other['mate_ids'] and row['partner']['chrom']==other['chrom'] and row['partner']['pos']==other['pos'] and other['partner']['chrom']==row['chrom'] and other['partner']['pos']==row['pos'])
                    if reciprocal and row['partner']['local_side'] is not None and other['partner']['local_side'] is not None:
                        reciprocal=(row['partner']['local_side']==other['partner']['remote_side'] and row['partner']['remote_side']==other['partner']['local_side'])
                    if reciprocal:
                        row['mate_status']='validated_explicit_pair';row['source_records']+=other['source_records'];row['mate_evidence']=other
                        db.execute('UPDATE records SET consumed=1 WHERE ordinal=?',(other_ord,))
                    else:row['mate_status']='unresolved_mate'
                else:row['mate_status']='orphan' if not candidates else 'ambiguous_mates'
            counts['events']+=1;counts[row['mate_status']]+=1
            if row['gene_pair']:counts['annotated_gene_pairs']+=1
            writer.writerow([row['id'],caller,json.dumps(row['gene_pair']),row['mate_status'],json.dumps(row)])
            # Prefer interpretable pairs, then native order. No support-derived score.
            review.append(row);review.sort(key=lambda r:(r['gene_pair'] is None,r['gene_pair'] or [],r['chrom'],r['pos'],r['id']));del review[100:]
    db.close();(out/'.sv.sqlite').unlink()
    summary=dict(schema_version=1,sample=sample,build='hg38',caller=caller,status='completed',counts=dict(counts),records=review,total_records=counts['events'],truncated=counts['events']>len(review),review_limit=100,reference_audit=dict(audit,gff_source=str(gff),fai_source=str(fai)),clinical_annotation_status='unavailable',files=[p.name for p in sorted(out.glob(prefix+'*'))])
    (out/(prefix+'.sv_annotation.json')).write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('sample','caller','vcf','gff','fai','output'):p.add_argument('--'+name,required=True)
    annotate(**vars(p.parse_args()))
