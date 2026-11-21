#!/usr/bin/env python3
"""Exact two-caller, PASS-only consensus across exact target regions; never infer somatic origin."""
import argparse
from collections import Counter, defaultdict
from bisect import bisect_left
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote
import pysam
import pysam.bcftools
from check_reference import contigs, contig_aliases


def checksum(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def resolve_queries(config, gff, reference, fai):
    genes = json.loads(Path(config).read_text())['genes']
    lengths = contigs(fai)
    aliases = contig_aliases(reference, lengths)
    def chrom(name):
        found = name if name in lengths else aliases.get(name)
        if found is None: raise ValueError('Query contig absent from reference: ' + name)
        return found
    wanted = set(genes['snv'].get('pathogenic', [])) - set(genes['snv'].get('pharmacogenomics', []))
    # Three streaming passes avoid retaining a whole-genome GFF in memory.
    def features():
        with open(gff) as stream:
            for line in stream:
                if line.startswith('#'): continue
                f = line.rstrip().split('\t')
                if len(f) != 9: continue
                attrs = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
                attrs = {k: unquote(v) for k, v in attrs.items()}
                yield f, attrs
    ids = {}
    for f, a in features():
        name = a.get('Name') or a.get('gene')
        if f[2] == 'gene' and name in wanted:
            if name in ids: raise ValueError('Ambiguous gene annotation: ' + name)
            ids[name] = a['ID']
    transcripts = defaultdict(list)
    for f, a in features():
        if f[2] in ('mRNA', 'transcript'):
            for gene, gid in ids.items():
                if gid in a.get('Parent', '').split(','): transcripts[gene].append(a['ID'])
    candidates = {t for ts in transcripts.values() for t in ts}
    cds = defaultdict(list)
    for f, a in features():
        if f[2] == 'CDS':
            for parent in a.get('Parent', '').split(','):
                if parent in candidates: cds[parent].append((f[0], int(f[3]), int(f[4]), f[6]))
    queries = []
    for gene in sorted(wanted):
        cfg = genes['snv']['transcripts'].get(gene)
        if not cfg or not cfg.get('variants'): raise ValueError('Missing pathogenic SNV query positions: ' + gene)
        choices = transcripts[gene]
        t = cfg.get('transcript')
        if t is None and choices:
            # Fail on ties rather than reproduce NASVAR's potentially ambiguous hash iteration.
            sizes={x:sum(e-s+1 for _,s,e,_ in cds[x]) for x in choices}
            best=[x for x in choices if sizes[x]==max(sizes.values())]
            if len(best)!=1:raise ValueError('Ambiguous longest CDS; configure transcript explicitly: '+gene)
            t=best[0]
        if t not in choices or not cds[t]: raise ValueError('Unresolved transcript: ' + gene + ':' + str(t))
        exons = sorted(set(cds[t]), key=lambda x: x[1])
        names = {x[0] for x in exons}; strands = {x[3] for x in exons}
        if len(names) != 1 or len(strands) != 1 or next(iter(strands)) not in ('+', '-'):
            raise ValueError('Inconsistent CDS: ' + t)
        strand = exons[0][3]; c = chrom(exons[0][0])
        if any(s < 1 or e < s or e > lengths[c] for _,s,e,_ in exons): raise ValueError('Invalid CDS bounds: ' + t)
        if any(a[2] >= b[1] for a,b in zip(exons,exons[1:])): raise ValueError('Overlapping CDS: ' + t)
        for label, positions in sorted(cfg['variants'].items()):
            for pos in positions:
                remaining = int(pos); genomic = None
                if remaining < 1: raise ValueError('CDS positions must be 1-based')
                for _, start, end, _ in (exons if strand == '+' else reversed(exons)):
                    size = end-start+1
                    if remaining <= size:
                        genomic = start+remaining-1 if strand == '+' else end-remaining+1
                        break
                    remaining -= size
                if genomic is None: raise ValueError('Query position exceeds CDS: ' + gene + ':' + str(pos))
                queries.append(dict(kind='snv', gene=gene, label=label, transcript=t, cds_position=pos,
                                    chrom=c, start=genomic-1, end=genomic, strand=strand))
    for gene, cfg in sorted(genes.get('itd', {}).items()):
        c = chrom(cfg['chrom']); start, end = int(cfg['start'])-1, int(cfg['end'])
        if not 0 <= start < end <= lengths[c]: raise ValueError('Invalid indel window: ' + gene)
        queries.append(dict(kind='indel', gene=gene, label=cfg.get('label', 'ITD query window'),
                            transcript=cfg.get('transcript'), chrom=c, start=start, end=end,
                            nasvar_thresholds={k:v for k,v in cfg.items() if k.startswith('min_')}))
    if not queries: raise ValueError('Empty NASVAR consensus query set')
    return queries


def target_intervals(path, lengths):
    regions=defaultdict(list)
    with open(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith(('#','track ','browser ')):continue
            f=line.split();chrom=f[0];start,end=int(f[1]),int(f[2])
            if chrom not in lengths or not 0<=start<end<=lengths[chrom]:raise ValueError('Invalid target interval: '+line.strip())
            regions[chrom].append((start,end))
    merged={}
    for chrom,intervals in regions.items():
        out=[]
        for start,end in sorted(intervals):
            if out and start<=out[-1][1]:out[-1]=(out[-1][0],max(end,out[-1][1]))
            else:out.append((start,end))
        merged[chrom]=out
    if not merged:raise ValueError('Empty target BED')
    return merged


def in_targets(regions, chrom, start, end):
    intervals=regions.get(chrom,[])
    i=bisect_left(intervals,(end,))-1
    return i>=0 and intervals[i][1]>start


def evidence(rec, sample):
    values = rec.samples[sample]
    return dict(filter=list(rec.filter) or ['.'], qual=rec.qual,
                **{key: values.get(key) for key in ('DP', 'AD', 'AF', 'VAF', 'GT') if key in values})


def consensus(sample, sources, fasta, fai, config, gff, reference, output, image=None, targets=None, query_annotations=None):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',sample): raise ValueError('Unsafe sample identifier')
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    if set(sources) != {'deepsomatic','clairsto'}: raise ValueError('Both somatic callers are required')
    lengths=contigs(fai)
    if targets is None:raise ValueError('Exact target BED required for somatic consensus')
    regions=target_intervals(targets,lengths)
    queries=[];annotation_error=None
    if all([config,gff,reference]):
        try:queries=resolve_queries(config,gff,reference,fai)
        except ValueError as exc:annotation_error=str(exc)
    if query_annotations:
        queries=json.loads(Path(query_annotations).read_text())['queries']
        if any(q['chrom'] not in lengths or not 0<=q['start']<q['end']<=lengths[q['chrom']] for q in queries):raise ValueError('Invalid recorded query annotation')
    seen = defaultdict(lambda: defaultdict(list)); excluded = Counter(); normalized_paths=[]
    for caller, paths in sorted(sources.items()):
        for i, path in enumerate(paths):
            normalized = out / f'.{caller}.{i}.vcf'
            # -c e rejects reference mismatches; -m -any splits, but does not atomize MNVs.
            pysam.bcftools.norm('-f',str(fasta),'-c','e','-m','-any','-Ov','-o',str(normalized),str(path),catch_stdout=False)
            normalized_paths.append((caller,normalized))
            with pysam.VariantFile(str(normalized)) as vcf:
                if list(vcf.header.samples) != [sample]: raise ValueError('VCF sample mismatch: ' + str(path))
                for rec in vcf:
                    if len(rec.alts or ()) != 1: excluded[caller+':missing_allele'] += 1; continue
                    alt = rec.alts[0]; ref = rec.ref
                    if not re.fullmatch('[ACGT]+',ref) or not re.fullmatch('[ACGT]+',alt) or alt == ref:
                        excluded[caller+':unsupported_allele'] += 1; continue
                    pure_indel = ((len(ref)==1 and len(alt)>1 and (alt.startswith(ref) or alt.endswith(ref))) or
                                  (len(alt)==1 and len(ref)>1 and (ref.startswith(alt) or ref.endswith(alt))))
                    kind = 'snv' if len(ref)==len(alt)==1 else 'indel' if pure_indel else None
                    if kind is None: excluded[caller+':complex_substitution'] += 1; continue
                    matching = [q for q in queries if q['kind']==kind and q['chrom']==rec.contig and rec.start<q['end'] and rec.start+len(ref)>q['start']]
                    if not in_targets(regions,rec.contig,rec.start,rec.start+len(ref)):
                        excluded[caller+':outside_targets'] += 1; continue
                    key = (rec.contig,rec.pos,ref,alt)
                    ev = evidence(rec,sample)
                    # Retain only PASS candidate alleles and small optional NASVAR-query evidence.
                    if ev['filter']==['PASS'] or matching:
                        if ev not in seen[key][caller]:seen[key][caller].append(ev)
    # Collect non-PASS counterpart evidence without retaining every emitted reference record.
    for caller,path in normalized_paths:
        with pysam.VariantFile(str(path)) as vcf:
            for rec in vcf:
                if len(rec.alts or ())!=1:continue
                key=(rec.contig,rec.pos,rec.ref,rec.alts[0])
                if key in seen:
                    ev=evidence(rec,sample)
                    if ev not in seen[key][caller]:seen[key][caller].append(ev)
        path.unlink()
    lengths=contigs(fai); order={c:i for i,c in enumerate(lengths)}
    rows=[]; counts=Counter(snv=0,indel=0)
    headers=pysam.VariantHeader()
    headers.add_meta('source',value='FrankONTstein two-caller tumor-only somatic candidates')
    for c,n in lengths.items(): headers.contigs.add(c,length=n)
    headers.add_meta('INFO',items=[('ID','CALLERS'),('Number','.'),('Type','String'),('Description','Distinct exact-PASS supporting callers')])
    headers.add_meta('INFO',items=[('ID','GENES'),('Number','.'),('Type','String'),('Description','NASVAR query genes; not a pathogenicity assertion')])
    for caller, tag in [('deepsomatic','DS'),('clairsto','CS')]:
        for field,typ in [('DP','Integer'),('AF','Float')]:
            headers.add_meta('INFO',items=[('ID',tag+'_'+field),('Number',1),('Type',typ),('Description',caller+' reported '+field+'; missing if inconsistent across duplicates')])
    headers.formats.add('GT',1,'String','Consensus genotype is not inferred')
    headers.add_sample(sample)
    writers={k:pysam.VariantFile(str(out/f'{sample}.consensus.{k}.vcf.gz'),'wz',header=headers) for k in ('snv','indel')}
    for key, ev in sorted(seen.items(),key=lambda kv:(order[kv[0][0]],*kv[0][1:])):
        c,pos,ref,alt=key; kind='snv' if len(ref)==len(alt)==1 else 'indel'
        qs=[q for q in queries if q['kind']==kind and q['chrom']==c and pos-1<q['end'] and pos-1+len(ref)>q['start']]
        support=[caller for caller in sorted(sources) if any(e['filter']==['PASS'] for e in ev.get(caller,[]))]
        status='consensus' if len(support)==2 else 'single_caller_pass' if support else 'no_pass_support'
        rows.append(dict(chrom=c,pos=pos,ref=ref,alt=alt,kind=kind,genes=sorted({q['gene'] for q in qs}),
                         queries=qs,status=status,pass_callers=support,evidence=dict(ev)))
        if status!='consensus': continue
        counts[kind]+=1
        rec=writers[kind].new_record(contig=c,start=pos-1,alleles=(ref,alt),filter='PASS')
        rec.info['CALLERS']=tuple(support)
        if qs:rec.info['GENES']=tuple(sorted({q['gene'] for q in qs}))
        rec.samples[sample]['GT']=(None,None)
        for caller,tag in [('deepsomatic','DS'),('clairsto','CS')]:
            for field in ('DP','AF'):
                vals=[]
                for e in ev[caller]:
                    if e['filter']!=['PASS']:continue
                    value=e.get(field,e.get('VAF') if field=='AF' else None)
                    if isinstance(value,(tuple,list)): value=value[0] if len(value)==1 else None
                    if value is not None: vals.append(value)
                if vals and len(set(vals))==1:rec.info[tag+'_'+field]=vals[0]
        writers[kind].write(rec)
    for kind, writer in writers.items():
        writer.close();pysam.tabix_index(str(out/f'{sample}.consensus.{kind}.vcf.gz'),preset='vcf',force=True)
    with (out/f'{sample}.consensus.evidence.tsv').open('w') as f:
        keys=['chrom','pos','ref','alt','kind','genes','status','pass_callers','queries','evidence']
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()})
    agreement=Counter(r['status'] for r in rows)
    caller_pass={caller:{kind:sum(r['kind']==kind and caller in r['pass_callers'] for r in rows) for kind in ('snv','indel')} for caller in sources}
    review=[r for r in rows if r['status']=='consensus'][:100]
    review += [r for r in rows if r['status']!='consensus'][:100]
    for r in rows:
        if r['queries'] and r not in review:review.append(r)
    summary=dict(schema_version=2,sample=sample,status='completed',scope='targets',counts=dict(counts),queries=queries,
                 annotation_status='unavailable' if annotation_error or not queries else 'available',annotation_error=annotation_error,
                 evidence=review,evidence_records=len(rows),evidence_truncated=len(review)<len(rows),review_limit_per_status=100,
                 caller_pass_counts=caller_pass,agreement_counts=dict(agreement),excluded_records=dict(excluded),
                 target_bases=sum(end-start for intervals in regions.values() for start,end in intervals),
                 policy='Exact normalized allele; exact PASS from both DeepSomatic and ClairS-TO across targets; tumor-only candidate, not proven somatic',
                 boundary_policy='Exact target BED, no padding; normalized REF-span overlap including anchor; 0-based half-open intervals',
                 nasvar_threshold_policy='NASVAR queries are optional annotations only; no NASVAR gate or thresholds on consensus; PGx interpretation remains separate',
                 query_status=[dict(query=q,status='consensus' if any(r['status']=='consensus' and q in r['queries'] for r in rows) else 'no_consensus_detected') for q in queries])
    (out/f'{sample}.consensus.json').write_text(json.dumps(summary,indent=2)+'\n')
    provenance=dict(image=image,pysam=pysam.__version__,htslib=pysam.__samtools_version__,
                    checksums={str(p):checksum(p) for p in [p for p in [targets,query_annotations,config,gff,reference,fai,fasta] if p is not None]},
                    fasta=str(fasta),sources={c:[dict(path=str(p),sha256=checksum(p)) for p in ps] for c,ps in sources.items()})
    (out/f'{sample}.consensus.provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('sample','fasta','fai','targets','output'):p.add_argument('--'+name,required=True)
    for name in ('config','gff','reference'):p.add_argument('--'+name)
    p.add_argument('--query-annotations',help='Recorded query JSON used only for optional locus annotations')
    p.add_argument('--deepsomatic',nargs='+',required=True);p.add_argument('--clairsto',nargs='+',required=True);p.add_argument('--image')
    a=p.parse_args();consensus(a.sample,dict(deepsomatic=a.deepsomatic,clairsto=a.clairsto),a.fasta,a.fai,a.config,a.gff,a.reference,a.output,a.image,a.targets,a.query_annotations)
