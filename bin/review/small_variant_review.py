#!/usr/bin/env python3
"""Precompute consolidated exact-allele review from complete annotation exports.

Selection uses the existing independent caller/callset/group budgets. Matching
caller evidence includes failed records; FILTER, DP, AD and AF remain native.
All selected alleles and their record evidence enter the report JSON.
"""
import argparse
from collections import defaultdict
import copy
import csv
import gzip
import json
from pathlib import Path
from report_prioritization import small_candidate, consensus_candidate, sample_evidence, bound_candidates, POLICY_VERSION
from transcript_selection import load_panel, panel_provenance


def entries(path):
    with gzip.open(path,'rt') as stream:
        yield from csv.DictReader(stream,delimiter='\t')


def build(sample, paths, output, panel=None):
    csv.field_size_limit(100_000_000)
    panel_info=panel_provenance(panel);panel=load_panel(panel)
    grouped=defaultdict(list);sources=[];counts=[]
    for path in sorted(map(Path,paths)):
        metadata_path=path.with_name(path.name.replace('.annotation.tsv.gz','.annotation.json'))
        meta=json.loads(metadata_path.read_text())
        if meta.get('sample')!=sample or meta.get('build')!='hg38':raise ValueError('Small review sample/build mismatch')
        sources.append((path,meta));ordinary=[];broad=[];total=0;broad_total=0
        for entry in entries(path):
            if entry['filter']!='PASS':continue
            ann=json.loads(entry['annotation']);allele=json.loads(entry['allele'])
            candidate=small_candidate(ann,['PASS'],panel)
            if candidate:
                total+=1
                reasons=[candidate['reason']]
                if candidate['reason']=='ClinVar pathogenic / likely pathogenic':
                    secondary=copy.deepcopy(ann);secondary['clinvar']=[]
                    consequence=small_candidate(secondary,['PASS'],panel)
                    if consequence:
                        if consequence['reason'] not in reasons:reasons.append(consequence['reason'])
                        candidate['trigger_transcripts']=candidate.get('trigger_transcripts',[])+[t for t in consequence.get('trigger_transcripts',[]) if t not in candidate.get('trigger_transcripts',[])]
                ordinary.append(dict(candidate,coordinates=allele,caller=meta['caller'],kind=meta['kind'],inclusion_reasons=reasons))
            if meta['caller']=='consensus':
                candidate=consensus_candidate(ann,['PASS'],panel)
                if candidate:
                    broad_total+=1
                    if True:
                        tx=candidate['display'].get('transcript') or {}
                        broad.append(dict(candidate,coordinates=allele,caller='consensus',kind=meta['kind'],priority=1,review_group='consensus_review',
                                          gene=tx.get('SYMBOL'),gene_context=candidate['display'],transcript=tx.get('Feature'),consequence=tx.get('Consequence'),
                                          hgvsc=tx.get('HGVSc'),hgvsp=tx.get('HGVSp'),inclusion_reasons=['Consensus review: '+candidate['reason']]))
        ordinary=bound_candidates(ordinary)
        for row in ordinary+broad:grouped[tuple(row['coordinates'])].append(row)
        counts.append(dict(caller=meta['caller'],kind=meta['kind'],selected_total=total,selected_display=len(ordinary),consensus_review_total=broad_total,consensus_review_display=len(broad)))
    evidence=defaultdict(list)
    for path,meta in sources:
        for entry in entries(path):
            allele=tuple(json.loads(entry['allele']))
            if allele not in grouped:continue
            ann=json.loads(entry['annotation'])
            row=dict(caller=meta['caller'],kind=meta['kind'],filter=entry['filter'],
                     sample_evidence=sample_evidence(entry.get('source_record',''),meta['caller']),source_record=entry.get('source_record'),
                     clinvar=ann.get('clinvar',[]),population=ann.get('population',{}),cosmic=ann.get('cosmic',{}),native_records=json.loads(entry.get('native_records','[]')),file=path.name)
            if row not in evidence[allele]:evidence[allele].append(row)
    rows=[dict(coordinates=list(key),selections=selections,evidence=evidence[key],
               exact_pass_agreement={'clairsto','deepsomatic'} <= {e['caller'] for e in evidence[key] if e['filter']=='PASS'})
          for key,selections in sorted(grouped.items(),key=lambda kv:(min(r['priority'] for r in kv[1]),kv[1][0].get('gene') or '',kv[0]))]
    result=dict(schema_version=1,sample=sample,build='hg38',status='completed',report_policy_version=POLICY_VERSION,report_panel=panel_info,
                records=rows,total_records=len(rows),selection_counts=counts,files=[p.name for p,_ in sources],
                policy='COSMIC blood/lymphoid primary review; other COSMIC and ClinVar/VEP candidates in expanded review. Caller PASS and population eligibility required. Exact normalized allele grouping; native FILTER/DP/AD/AF retained; no caller-count or AF pooling.')
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',required=True);p.add_argument('--inputs',nargs='+',required=True);p.add_argument('--output',required=True);p.add_argument('--panel')
    a=p.parse_args();build(a.sample,a.inputs,a.output,a.panel)
