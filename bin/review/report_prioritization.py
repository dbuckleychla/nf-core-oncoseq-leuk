"""Scientific selection for concise reports, applied before bounded rendering.

These are research review criteria, not ACMG classification or somatic inference.
Full annotation exports are never filtered by this module.
"""
import csv
import math
import re

from transcript_selection import load_panel, select_transcript, effects, matches_preference

POLICY_VERSION = 3
LOSS_OF_FUNCTION = {'transcript_ablation', 'frameshift_variant', 'stop_gained',
                    'splice_acceptor_variant', 'splice_donor_variant'}


def text_values(value):
    if isinstance(value, (list, tuple)):
        return '|'.join(map(str, value))
    return str(value or '')


def sample_evidence(source_record, caller):
    """Extract per-caller depth/AF; calculate AD fraction only when native AF is absent."""
    fields=source_record.strip().split('\t')
    if len(fields)<8:return {}
    def number(value,integer=False):
        if value in (None,'','.'):return None
        if ',' in str(value):return None
        try:n=int(value) if integer else float(value)
        except (TypeError,ValueError):return None
        return n if math.isfinite(n) and n>=0 else None
    if caller=='consensus':
        info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
        return {name:dict(depth=number(info.get(tag+'_DP'),True),af=number(info.get(tag+'_AF')),af_source='consensus native caller INFO')
                for name,tag in [('deepsomatic','DS'),('clairsto','CS')]}
    fmt=dict(zip(fields[8].split(':'),fields[9].split(':'))) if len(fields)>9 else {}
    af=None;source=None
    for tag in ('AF','VAF'):
        af=number(fmt.get(tag))
        if af is not None:
            source='FORMAT/'+tag
            break
    info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
    split_multiallelic=',' in info.get('ONCOSEQ_REVIEW_ORIG','')
    if af is None and not split_multiallelic:
        ad=[number(x,True) for x in fmt.get('AD','').split(',')]
        if len(ad)==2 and all(x is not None for x in ad) and sum(ad)>0:
            af=ad[1]/sum(ad);source='FORMAT/AD alt / (ref + alt)'
    if af is not None and af>1:raise ValueError('Sample AF outside [0,1]')
    ad=fmt.get('AD','').split(',')
    support=number(ad[1],True) if ',' not in fields[4] and len(ad)==2 else None
    return {caller:dict(depth=number(fmt.get('DP'),True),af=af,af_source=source,
                        supporting_reads=support,support_source='FORMAT/AD ALT count' if support is not None else None)}


def small_candidate(annotation, filters, panel=None):
    panel = load_panel() if panel is None else panel
    if filters != ['PASS'] or annotation['population']['status'] == 'excluded_common':
        return None
    assertions = annotation.get('clinvar', [])
    sig = '|'.join(text_values(c.get('CLNSIG')) for c in assertions).lower()
    conflict = 'conflict' in sig or any(c.get('CLNSIGCONF') or 'conflicting' in text_values(c.get('CLNREVSTAT')).lower() for c in assertions)
    benign = 'benign' in sig
    uncertain = 'uncertain' in sig
    reviewed_pathogenic = any(
        re.search(r'(^|[|/,])(?:likely_)?pathogenic($|[|/,])', text_values(c.get('CLNSIG')).lower())
        and re.search(r'(^|[|,])(?:criteria_provided|reviewed_by_expert_panel|practice_guideline)(?:$|[|,_])', text_values(c.get('CLNREVSTAT')).lower())
        for c in assertions)
    pathogenic = bool(reviewed_pathogenic and not (conflict or benign or uncertain))
    # VEP HIGH on an incidental/noncoding transcript is insufficient.
    transcripts = [t for t in annotation.get('consequences', [])
                   if t.get('SYMBOL') and t.get('BIOTYPE') == 'protein_coding'
                   and (t.get('MANE_SELECT') or t.get('MANE_PLUS_CLINICAL') or t.get('CANONICAL') == 'YES'
                        or matches_preference(t, panel.get(t.get('SYMBOL'), '')))
                   and not t.get('FLAGS')]
    damaging = [t for t in transcripts if LOSS_OF_FUNCTION.intersection(t.get('Consequence', '').split('&'))
                and t.get('IMPACT') == 'HIGH']
    coding = [t for t in transcripts if t['SYMBOL'] in panel and effects(t).intersection({'missense_variant', 'inframe_insertion', 'inframe_deletion', 'protein_altering_variant'})]
    strong = bool(damaging and not (conflict or benign))
    if not pathogenic and not strong and not coding:
        return None
    choices = transcripts if pathogenic else damaging if strong else coding
    if not choices:  # Keep reviewed noncoding pathogenic variants with a named gene.
        choices = [t for t in annotation.get('consequences', []) if t.get('SYMBOL')]
    if not choices:
        return None
    selection = select_transcript(choices, panel)
    tx = selection['transcript']
    group = 'supporting' if benign or tx['SYMBOL'] not in panel else 'findings' if pathogenic or strong else 'coding'
    reason = 'ClinVar pathogenic / likely pathogenic' if pathogenic else 'Predicted loss of function' if strong else 'Panel coding variant; significance not established'
    return dict(reason=reason, review_group=group, panel_gene=tx['SYMBOL'] in panel,
                clinvar=assertions,
                gene_context=selection, trigger_transcripts=damaging if strong and not pathogenic else [],
                priority=0 if pathogenic else 1 if strong else 2, gene=tx['SYMBOL'], transcript=tx.get('Feature'),
                consequence=tx.get('Consequence'), hgvsc=tx.get('HGVSc'), hgvsp=tx.get('HGVSp'),
                population_af=annotation['population'].get('combined_af'),
                clinvar_significance=sig or None,
                clinvar_review='; '.join(sorted({text_values(c.get('CLNREVSTAT')) for c in assertions})))


def consensus_candidate(annotation, filters, panel=None):
    """Broad review of exact PASS consensus: any HIGH transcript or P/LP assertion."""
    if filters != ['PASS']:return None
    assertions=annotation.get('clinvar',[])
    pathogenic=any(re.search(r'(^|[|/,])(?:likely_)?pathogenic($|[|/,])',text_values(c.get('CLNSIG')).lower()) for c in assertions)
    high=[t for t in annotation.get('consequences',[]) if t.get('IMPACT')=='HIGH']
    if not high and not pathogenic:return None
    txs=high or annotation.get('consequences',[])
    return dict(reason=' + '.join(x for x,yes in [('VEP HIGH',bool(high)),('ClinVar P/LP',pathogenic)] if yes),
        display=select_transcript(annotation.get('consequences', []), panel),
        transcripts=[{k:t.get(k) for k in ['SYMBOL','Feature','MANE_SELECT','HGVSc','HGVSp','Consequence','IMPACT']} for t in txs],
        clinvar=assertions,population=annotation.get('population',{}))


def candidate_sort(row):
    return (row['priority'], row['gene'], row['coordinates'])


def bound_candidates(rows, limit=None):
    # Independent budgets keep additional coding review from displacing findings.
    return [r for group in ('findings', 'coding', 'supporting')
            for r in sorted((r for r in rows if r.get('review_group', 'findings') == group),
                            key=candidate_sort)[:limit]]


def fusion_pairs(path):
    with open(path) as stream:
        return {tuple(sorted((r['gene'], r['partner']))) for r in csv.DictReader(stream, delimiter='\t')
                if r['gene'] and r['partner'] and r['gene'] != r['partner']}


def fusion_candidate(row, known_pairs):
    """Resolve curated pairs across overlapping genes without changing native labels."""
    if row['caller'] == 'NASVAR':
        e = row['evidence']
        support = e.get('supporting_reads')
        a={g['name'] for g in row.get('genes1',[])}
        b={g['name'] for g in row.get('genes2',[])}
        pairs={tuple(sorted((x,y))) for x in a for y in b if x!=y}
        if row.get('gene_pair'):pairs.add(tuple(sorted(row['gene_pair'])))
        svtype='native fusion'
    else:
        if row.get('filter') != ['PASS']:
            return None
        # Both explicit mates must pass; do not promote a failed counterpart.
        if row.get('mate_evidence') and row['mate_evidence'].get('filter') != ['PASS']:
            return None
        e = row['evidence']; info = e.get('info', {})
        svtype = info.get('SVTYPE')
        if svtype and svtype not in ('BND', 'TRA', 'INV', 'DEL', 'DUP', 'DUP:TANDEM'):
            return None
        if svtype in ('INV','DEL','DUP','DUP:TANDEM'):
            partner=row.get('partner') or {}
            if partner.get('chrom')!=row.get('chrom') or not isinstance(partner.get('pos'),int) or not isinstance(row.get('pos'),int) or not 1<=row['pos']<partner['pos']:
                return None
        if not svtype and row.get('partner', {}).get('encoding') not in ('BND', 'CHR2/POS2'):
            return None
        a = {g['name'] for g in row.get('genes1', [])}
        b = {g['name'] for g in row.get('genes2', [])}
        pairs = {tuple(sorted((x, y))) for x in a for y in b if x != y}
        if row.get('gene_pair'):pairs.add(tuple(sorted(row['gene_pair'])))
        fmt = e.get('format', {}).get(row.get('vcf_sample',row['sample']), {})
        support = fmt.get('DV')
        if support is None:support = info.get('SUPPORT', info.get('SR'))
        mate = row.get('mate_evidence')
        if mate:
            mi = mate['evidence'].get('info', {})
            ms = mate['evidence'].get('format', {}).get(row.get('vcf_sample',row['sample']), {}).get('DV')
            if ms is None:ms = mi.get('SUPPORT', mi.get('SR'))
            if not isinstance(ms, (float, int)) or not math.isfinite(ms):return None
            if isinstance(support, (float, int)):support = min(support, ms)
    matches = pairs & known_pairs
    if len(matches) != 1 or not isinstance(support, (float, int)) or not math.isfinite(support) or support < 3:
        return None
    structural=svtype in ('INV','DEL','DUP','DUP:TANDEM')
    return dict(pair=list(next(iter(matches))), caller=row['caller'], support=support,
                native_id=row.get('native_id'), chrom=row['chrom'], pos=row['pos'],
                partner=row['partner'], coordinate_convention=row['coordinate_convention'],
                overlapping_gene_ambiguity=len(row.get('genes1',[]))>1 or len(row.get('genes2',[]))>1 or not bool(row.get('gene_pair')),
                reason=('Curated gene-pair structural evidence' if structural else 'Curated fusion pair')+'; at least 3 supporting reads',
                evidence_kind='supporting_structural' if structural else 'fusion_candidate',native_type=svtype,
                genes1=row.get('genes1',[]),genes2=row.get('genes2',[]),native_gene_pair=row.get('native_gene_pair') if row['caller']=='NASVAR' else [info['GENE1'],info['GENE2']] if info.get('GENE1') and info.get('GENE2') else None,
                evidence=row['evidence'],filter=row.get('filter'),source_records=row.get('source_records',[]),mate_evidence=row.get('mate_evidence'),
                breakpoint_annotations=row.get('breakpoint_annotations',[]),mate_status=row.get('mate_status'))
