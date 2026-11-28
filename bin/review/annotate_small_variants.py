#!/usr/bin/env python3
"""Disk-backed exact-allele annotation; caller FILTER and evidence are preserved."""
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
from somatic_consensus import target_intervals, in_targets
from check_reference import contigs
from indexed_reference import IndexedReference
from progress import progress


def key(rec):
    return json.dumps([rec.contig, rec.pos, rec.ref, rec.alts[0]], separators=(',', ':'))


def prepare(sources, fasta, targets, output):
    out=Path(output);out.mkdir(parents=True, exist_ok=True)
    db=sqlite3.connect(out/'catalog.sqlite')
    db.executescript('CREATE TABLE alleles (key TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, candidate INTEGER, annotation TEXT); CREATE TABLE records (source INTEGER, key TEXT, record TEXT); CREATE INDEX record_key ON records(key); CREATE INDEX record_source ON records(source,key); CREATE INDEX allele_chrom ON alleles(chrom,pos);')
    db.executescript('CREATE TABLE original_records (source INTEGER, site TEXT, record TEXT); CREATE INDEX original_site ON original_records(source,site); CREATE TABLE native_records (source INTEGER, key TEXT, id TEXT, native_json TEXT); CREATE INDEX native_key ON native_records(source,key,id);')
    lengths=contigs(str(fasta)+'.fai');regions=target_intervals(targets,lengths)
    rows=json.loads(Path(sources).read_text());excluded=Counter()
    for i, source in enumerate(rows):
        source['target_filter_counts']=dict(pass_records=0,unassessed_records=0,failed_records=0)
        source['excluded_records']=dict(unsupported_allele=0,outside_targets=0)
        with open(source['path'],'rb') as native_stream,pysam.VariantFile(native_stream) as native:
            for original in native:
                site=json.dumps([original.contig,original.pos,original.ref,list(original.alts or [])],separators=(',',':'))
                db.execute('INSERT INTO original_records VALUES (?,?,?)',(i,site,str(original)))
        db.commit()
        normalized=out/f'{i}.normalized.vcf'
        pysam.bcftools.norm('-f',str(fasta),'-c','e','-m','-any','--old-rec-tag','ONCOSEQ_REVIEW_ORIG','-Ov','-o',str(normalized),source['path'],catch_stdout=False)
        with pysam.VariantFile(normalized) as vcf:
            if list(vcf.header.samples)!=[source['sample']]:raise ValueError('Caller sample mismatch: '+source['path'])
            (out/f'{i}.header.vcf').write_text(str(vcf.header))
            for rec in vcf:
                if len(rec.alts or ())!=1 or not re.fullmatch('[ACGT]+',rec.ref) or not re.fullmatch('[ACGT]+',rec.alts[0]) or rec.ref==rec.alts[0]:
                    excluded['unsupported_allele']+=1;source['excluded_records']['unsupported_allele']+=1;continue
                if not in_targets(regions,rec.contig,rec.start,rec.start+len(rec.ref)):
                    excluded['outside_targets']+=1;source['excluded_records']['outside_targets']+=1;continue
                filter_key='pass_records' if list(rec.filter)==['PASS'] else 'unassessed_records' if not list(rec.filter) else 'failed_records'
                source['target_filter_counts'][filter_key]+=1
                candidate=1
                k=key(rec)
                db.execute('INSERT INTO alleles VALUES (?,?,?,?,?,?,NULL) ON CONFLICT(key) DO UPDATE SET candidate=MAX(candidate,excluded.candidate)',(k,rec.contig,rec.pos,rec.ref,rec.alts[0],candidate))
                db.execute('INSERT INTO records VALUES (?,?,?)',(i,k,str(rec)))
                origin=rec.info.get('ONCOSEQ_REVIEW_ORIG')
                if isinstance(origin,tuple):origin=','.join(origin)
                if origin:
                    parts=origin.split('|')
                    site=json.dumps([parts[0],int(parts[1]),parts[2],parts[3].split(',')],separators=(',',':'))
                else:
                    site=json.dumps([rec.contig,rec.pos,rec.ref,list(rec.alts)],separators=(',',':'))
                originals=[text for (text,) in db.execute('SELECT record FROM original_records WHERE source=? AND site=?',(i,site))]
                matching=[text for text in originals if text.split('\t')[2]==(rec.id or '.')]
                db.execute('INSERT INTO native_records VALUES (?,?,?,?)',(i,k,rec.id,json.dumps(matching or originals)))
        normalized.unlink()
        db.commit()
    header=pysam.VariantHeader()
    for c,n in lengths.items():header.contigs.add(c,length=n)
    with pysam.VariantFile(out/'catalog.vcf','w',header=header) as writer:
        for c in lengths:
            for k,pos,ref,alt in db.execute('SELECT key,pos,ref,alt FROM alleles WHERE candidate=1 AND chrom=? ORDER BY pos,ref,alt',(c,)):
                writer.write(writer.new_record(contig=c,start=pos-1,alleles=(ref,alt),id=hashlib.sha256(k.encode()).hexdigest()[:24]))
    (out/'sources.json').write_text(json.dumps(rows))
    (out/'excluded.json').write_text(json.dumps(dict(excluded)))
    (out/'summary.json').write_text(json.dumps({'candidate_alleles':db.execute('SELECT COUNT(*) FROM alleles WHERE candidate=1').fetchone()[0]}))
    db.close()


def exact_records(vcf, rec):
    if isinstance(vcf, IndexedReference):return vcf.exact(rec)
    if rec.contig not in vcf.header.contigs:return []
    return [r for r in vcf.fetch(rec.contig,rec.start,rec.start+1) if r.pos==rec.pos and r.ref==rec.ref and r.alts==rec.alts]


def validate_annotation_contigs(handles, query_contigs):
    """Reject recognizable naming mismatches; never rename coordinates implicitly.

    Sparse/empty resource contigs may legitimately have no matching variants.
    Only fail when a queried chromosome is present under a different name.
    """
    def chromosome(name):
        name=name.removeprefix('chr')
        accession=re.fullmatch(r'NC_0*(\d+)\.\d+',name)
        if accession:
            number=int(accession[1])
            name={23:'X',24:'Y',12920:'MT'}.get(number,str(number))
        if name in ('M','MT'):return 'MT'
        return name if name in {str(i) for i in range(1,23)} | {'X','Y'} else None
    for resource,vcf in handles.items():
        declared=set(vcf.header.contigs)
        aliases={chromosome(c):c for c in declared if chromosome(c)}
        for contig in query_contigs:
            other=aliases.get(chromosome(contig))
            if contig not in declared and other:
                raise ValueError(f'Annotation chromosome naming mismatch in {resource}: query {contig}, resource {other}. Prepare matching reference contig names; refusing silent missing annotations.')


def counts(records):
    values=set()
    for rec in records:
        ac=rec.info.get('AC');an=rec.info.get('AN')
        if isinstance(ac,tuple):ac=ac[0] if len(ac)==1 else None
        if isinstance(an,tuple):an=an[0] if len(an)==1 else None
        if ac is None or an is None or an==0:continue
        if not isinstance(ac,int) or not isinstance(an,int) or not 0<=ac<=an:raise ValueError('Invalid population AC/AN')
        values.add((ac,an))
    if len(values)>1:raise ValueError('Conflicting population records for exact allele')
    if not values:return dict(ac=None,an=None,af=None)
    ac,an=next(iter(values));return dict(ac=ac,an=an,af=ac/an)


def population(genome,exome,rescue):
    usable=[d for d in (genome,exome) if d['an']]
    af=sum(d['ac'] for d in usable)/sum(d['an'] for d in usable) if usable else None
    status='retained_whitelist_rescue' if rescue else 'retained_frequency_unknown' if af is None else 'excluded_common' if af>0.01 else 'retained_rare'
    return dict(genome=genome,exome=exome,combined_af=af,status=status,whitelist_rescue=rescue)


def consequences(rec, fields):
    values=rec.info.get('CSQ',())
    return [dict(zip(fields,row.split('|'))) for row in values]


from transcript_selection import display_transcript, load_panel, panel_provenance


def review_group(annotation):
    if annotation['cosmic']['blood_match']:return 0
    if annotation['cosmic']['status'] == 'matched':return 1
    sig='|'.join(str(x.get('CLNSIG','')) for x in annotation['clinvar']).lower()
    conflict='conflict' in sig or ('pathogenic' in sig and ('benign' in sig or 'uncertain_significance' in sig)) or any('CLNSIGCONF' in x for x in annotation['clinvar'])
    pathogenic='pathogenic' in sig and not conflict
    benign='benign' in sig and not conflict
    if annotation['population']['whitelist_rescue'] or pathogenic:return 2
    if benign:return 5
    if any(x.get('IMPACT') in ['HIGH','MODERATE'] for x in annotation['consequences']):return 3
    return 4


def finish(catalog, vep, resources, sample, output, panel=None):
    panel_path=panel
    panel_info = panel_provenance(panel)
    panel = load_panel(panel)
    from report_prioritization import small_candidate, consensus_candidate, sample_evidence, POLICY_VERSION, bound_candidates
    cat=Path(catalog);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    db=sqlite3.connect(cat/'catalog.sqlite')
    handles={k:IndexedReference(v) for k,v in resources.items() if k in ['gnomad_genome','gnomad_exome','clinvar','clinvar_common_pathogenic'] or k == 'cosmic_targeted'}
    from cosmic_evidence import CosmicEvidence
    cosmic_evidence = CosmicEvidence(resources['cosmic_evidence'])
    for name,h in handles.items():
        required=['GENE','TRANSCRIPT'] if name.startswith('cosmic') else ['AC','AN'] if name.startswith('gnomad') else ['clinvar_cmn_path'] if name=='clinvar_common_pathogenic' else ['CLNSIG','CLNREVSTAT']
        if any(k not in h.header.info for k in required):raise ValueError('Missing annotation fields in '+name)
    validate_annotation_contigs(handles,[r[0] for r in db.execute('SELECT DISTINCT chrom FROM alleles WHERE candidate=1')])
    expected = db.execute('SELECT COUNT(*) FROM alleles WHERE candidate=1').fetchone()[0]
    processed = 0
    progress(f'Reference annotation started: {expected} unique alleles')
    with open(vep,'rb') as vep_stream,pysam.VariantFile(vep_stream) as vf:
        if 'CSQ' not in vf.header.info:raise ValueError('VEP did not produce CSQ header')
        fields=vf.header.info['CSQ'].description.split('Format: ')[-1].split('|')
        csq_header=vf.header.info['CSQ'].record
        for rec in vf:
            k=key(rec)
            if not db.execute('SELECT 1 FROM alleles WHERE key=? AND candidate=1',(k,)).fetchone():raise ValueError('Unexpected VEP allele')
            rescue=any(r.info.get('clinvar_cmn_path') in (1,(1,),True) for r in exact_records(handles['clinvar_common_pathogenic'],rec))
            clinical=[dict(record_id=r.id,**{k:r.info[k] for k in r.info if k.startswith('CLN') or k=='ALLELEID'}) for r in exact_records(handles['clinvar'],rec)]
            csq=consequences(rec,fields)
            ann=dict(population=population(counts(exact_records(handles['gnomad_genome'],rec)),counts(exact_records(handles['gnomad_exome'],rec)),rescue),clinvar=clinical,consequences=csq,display_transcript=display_transcript(csq, panel),csq=list(rec.info.get('CSQ',())))
            from cosmic import annotate as cosmic_annotate
            ann['cosmic']=cosmic_annotate(handles,rec,cosmic_evidence)
            ann['review_group']=review_group(ann)
            db.execute('UPDATE alleles SET annotation=? WHERE key=?',(json.dumps(ann),k))
            processed += 1
            if processed % 10000 == 0:progress(f'Reference annotation: {processed}/{expected}')
    db.commit()
    if db.execute('SELECT COUNT(*) FROM alleles WHERE candidate=1 AND annotation IS NULL').fetchone()[0]:raise ValueError('VEP omitted candidate alleles')
    progress(f'Reference annotation finished: {processed}/{expected}')
    for name, handle in handles.items():progress('Cache ' + name + ': ' + json.dumps(handle.statistics(), sort_keys=True))
    cosmic_evidence.close()
    sources=json.loads((cat/'sources.json').read_text())
    for i,source in enumerate(sources):
        caller=source['caller'];kind=source['kind'];dest=out/caller;dest.mkdir(exist_ok=True)
        progress('Writing annotated exports: ' + caller)
        prefix=f'{sample}.{caller}.{kind}.targets'
        header=pysam.VariantFile(cat/f'{i}.header.vcf').header.copy()
        header.add_record(csq_header)
        for name,num,typ,description in [('ONCOSEQ_REVIEW_POP_AF',1,'Float','Pooled valid gnomAD exome/genome AC/AN'),('ONCOSEQ_REVIEW_POP_STATUS',1,'String','Population eligibility; caller FILTER is unchanged'),('ONCOSEQ_REVIEW_CLINVAR','.', 'String','ClinVar assertions; percent-encoded JSON preserves conflicts')]:
            header.info.add(name,num,typ,description)
        header.info.add('ONCOSEQ_REVIEW_COSMIC', '.', 'String', 'Exact normalized COSMIC targeted COSV identifiers')
        header.info.add('ONCOSEQ_REVIEW_COSMIC_BLOOD', 0, 'Flag', 'Positive COSMIC targeted evidence in blood or lymphoid tissue; not sample somatic status')
        stats=Counter();review=[];total=0;report_candidates=[];report_total=0;report_groups=Counter();consensus_review=[];consensus_review_total=0
        contig_order={c:i for i,c in enumerate(header.contigs)}
        annotated=dest/(prefix+'.annotated.vcf.gz');prioritized=dest/(prefix+'.prioritized.vcf.gz')
        # Sorted by the original header contig order, then normalized position.
        with pysam.VariantFile(annotated,'wz',header=header) as aw,pysam.VariantFile(prioritized,'wz',header=header) as pw,gzip.open(dest/(prefix+'.annotation.tsv.gz'),'wt') as ts:
            writer=csv.writer(ts,delimiter='\t');writer.writerow(['allele','caller','kind','filter','population_status','annotation','source_record','native_records'])
            selected=cat/f'{i}.selected.vcf'
            with selected.open('w') as stream:
                stream.write(str(header))
                for chrom in header.contigs:
                    query='SELECT r.record FROM records r JOIN alleles a ON a.key=r.key WHERE r.source=? AND a.candidate=1 AND a.chrom=? ORDER BY a.pos,a.ref,a.alt,r.record'
                    for (text,) in db.execute(query,(i,chrom)):stream.write(text)
            with pysam.VariantFile(selected) as selected_vcf:
                for rec in selected_vcf:
                    k=key(rec);text=str(rec)
                    serialized=db.execute('SELECT annotation FROM alleles WHERE key=?',(k,)).fetchone()[0]
                    ann=json.loads(serialized)
                    rec.translate(header)
                    rec.info['CSQ']=tuple(ann['csq'])
                    pop=ann['population'];rec.info['ONCOSEQ_REVIEW_POP_STATUS']=pop['status']
                    if pop['combined_af'] is not None:rec.info['ONCOSEQ_REVIEW_POP_AF']=pop['combined_af']
                    from urllib.parse import quote
                    if ann['clinvar']:rec.info['ONCOSEQ_REVIEW_CLINVAR']=tuple(quote(json.dumps(c,separators=(',',':')),safe='') for c in ann['clinvar'])
                    cosmic_ids = sorted({identifier for source in ann['cosmic']['sources'] for identifier in source['identifiers']})
                    if cosmic_ids:rec.info['ONCOSEQ_REVIEW_COSMIC'] = tuple(cosmic_ids)
                    if ann['cosmic']['blood_match']:rec.info['ONCOSEQ_REVIEW_COSMIC_BLOOD'] = True
                    aw.write(rec);is_pass=list(rec.filter)==['PASS'];stats['annotated_records']+=1
                    stats[pop['status']]+=1
                    if is_pass:stats['pass_before_population']+=1
                    if is_pass and pop['status']!='excluded_common':pw.write(rec);stats['pass_after_population']+=1
                    native=[]
                    for (encoded_native,) in db.execute('SELECT native_json FROM native_records WHERE source=? AND key=? AND id IS ?',(i,k,rec.id)):
                        for original in json.loads(encoded_native):
                            if original not in native:native.append(original)
                    writer.writerow([k,caller,kind,','.join(rec.filter) or '.',pop['status'],serialized,text.strip(),json.dumps(native)])
                    total+=1
                    if total % 10000 == 0:progress(f'Export {caller}: {total} records')
                    if caller=='consensus':
                        selected_review=consensus_candidate(ann,list(rec.filter),panel)
                        if selected_review:
                            consensus_review_total+=1
                            if True:consensus_review.append(dict(selected_review,coordinates=json.loads(k),sample_evidence=sample_evidence(text,caller)))
                    candidate=small_candidate(ann,list(rec.filter),panel)
                    if candidate:
                        report_total+=1
                        report_groups[candidate['review_group']]+=1
                        report_candidates.append(dict(candidate,coordinates=json.loads(k),caller=caller,kind=kind,sample_evidence=sample_evidence(text,caller)))
                    if (is_pass or not list(rec.filter)) and pop['status']!='excluded_common':
                        row=dict(schema_version=1,sample=sample,build='hg38',record_type='small_variant',id=hashlib.sha256(k.encode()).hexdigest()[:24],coordinates=json.loads(k),contig_order=contig_order[rec.contig],coordinate_convention='VCF 1-based POS; normalized REF/ALT',caller=caller,kind=kind,filter=list(rec.filter) or ['.'],annotation=ann,source_records=[text.strip()])
                        review.append(row);review.sort(key=lambda r:(r['annotation']['review_group'],r['contig_order'],r['coordinates']))
                        del review[100:]
        report_candidates=bound_candidates(report_candidates)
        filters=source.get('target_filter_counts',{})
        unassessed=filters.get('unassessed_records',0)>0 and not (filters.get('pass_records',0) or filters.get('failed_records',0))
        if unassessed:prioritized.unlink()
        for p in ([annotated] if unassessed else [annotated,prioritized]):pysam.tabix_index(str(p),preset='vcf',force=True)
        stats.update({'unique_candidate_alleles':db.execute('SELECT COUNT(DISTINCT r.key) FROM records r JOIN alleles a ON a.key=r.key WHERE r.source=? AND a.candidate=1',(i,)).fetchone()[0]})
        summary=dict(schema_version=1,sample=sample,build='hg38',caller=caller,kind=kind,status='completed',pass_filter_status='not_applicable' if unassessed else 'applicable',source_filter_counts=filters,excluded_records=source.get('excluded_records',{}),counts=dict(stats),report_policy_version=POLICY_VERSION,report_candidates=report_candidates,report_total=report_total,report_groups=dict(report_groups),report_panel=panel_info,records=review,total_records=total,truncated=total>len(review),review_limit=100,files=[p.name for p in sorted(dest.glob(prefix+'*'))])
        if caller=='consensus':summary.update(consensus_review=consensus_review,consensus_review_total=consensus_review_total)
        (dest/(prefix+'.annotation.json')).write_text(json.dumps(summary,indent=2)+'\n')
    for h in handles.values():h.close()
    db.close()
    progress('Consolidating small-variant review evidence')
    from small_variant_review import build
    build(sample,list(out.glob('*/*.annotation.tsv.gz')),out/'annotation'/(sample+'.small_variant_review.json'),panel_path)


if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for name in ['sources','fasta','targets','output']:a.add_argument('--'+name,required=True)
    a=sub.add_parser('finish')
    for name in ['catalog','vep','resources','sample','output']:a.add_argument('--'+name,required=True)
    a.add_argument('--panel')
    args=vars(p.parse_args());command=args.pop('command')
    if command=='finish':args['resources']=json.loads(Path(args['resources']).read_text())
    globals()[command](**args)
