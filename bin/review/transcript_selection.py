"""Gene-aware transcript selection; never change allele or caller evidence."""
import csv
import hashlib
from pathlib import Path

DEFAULT_PANEL = Path(__file__).resolve().parents[2] / 'assets/reporting/leukemia-v1.tsv'
CONSEQUENCE_ORDER = 'transcript_ablation splice_acceptor_variant splice_donor_variant stop_gained frameshift_variant stop_lost start_lost transcript_amplification feature_elongation feature_truncation inframe_insertion inframe_deletion missense_variant protein_altering_variant splice_donor_5th_base_variant splice_region_variant splice_donor_region_variant splice_polypyrimidine_tract_variant incomplete_terminal_codon_variant start_retained_variant stop_retained_variant synonymous_variant coding_sequence_variant mature_miRNA_variant 5_prime_UTR_variant 3_prime_UTR_variant non_coding_transcript_exon_variant intron_variant NMD_transcript_variant non_coding_transcript_variant upstream_gene_variant downstream_gene_variant TFBS_ablation TFBS_amplification TF_binding_site_variant regulatory_region_ablation regulatory_region_amplification regulatory_region_variant intergenic_variant'.split()
FLANKING = {'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant'}
CODING = set(CONSEQUENCE_ORDER[:24])


def load_panel(path=None):
    source = Path(path or DEFAULT_PANEL)
    with source.open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or 'gene' not in reader.fieldnames:
            raise ValueError('Report panel requires a gene column')
        genes = {}
        for row in reader:
            gene = (row.get('gene') or '').strip()
            if not gene or gene in genes:
                raise ValueError('Empty or duplicate panel gene: ' + gene)
            genes[gene] = (row.get('preferred_transcript') or '').strip()
    return genes


def panel_provenance(path=None):
    source = Path(path or DEFAULT_PANEL)
    return dict(name=source.name, sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                genes=load_panel(source))


def effects(tx):
    return set(tx.get('Consequence', '').split('&'))


def context_rank(tx):
    terms = effects(tx)
    if terms & CODING:
        return 0
    return 2 if terms and terms <= FLANKING else 1


def matches_preference(tx, preference):
    return bool(preference) and preference in (
        tx.get('Feature'), tx.get('MANE_SELECT'), tx.get('MANE_PLUS_CLINICAL'),
        tx.get('HGVSc', '').split(':')[0])


def transcript_rank(tx, preference=''):
    return (not matches_preference(tx, preference), not bool(tx.get('MANE_SELECT')),
            not bool(tx.get('MANE_PLUS_CLINICAL')),
            not (tx.get('CANONICAL') == 'YES' and tx.get('BIOTYPE') == 'protein_coding'),
            min((CONSEQUENCE_ORDER.index(c) if c in CONSEQUENCE_ORDER else len(CONSEQUENCE_ORDER)
                 for c in effects(tx)), default=len(CONSEQUENCE_ORDER)), tx.get('Feature', ''))


def select_transcript(rows, panel=None):
    """Keep direct gene ambiguity; preferences never turn a flank into a coding hit."""
    panel = load_panel() if panel is None else panel
    named = [t for t in rows if t.get('SYMBOL')]
    if not named:
        tx = min(rows, key=lambda t: (context_rank(t), transcript_rank(t))) if rows else None
        return dict(transcript=tx, genes=[], ambiguous=False,
                    warnings=['Gene symbol unavailable'] if tx else [])
    best_context = min(context_rank(t) for t in named)
    direct = [t for t in named if context_rank(t) == best_context]
    genes = sorted({t['SYMBOL'] for t in direct})
    panel_genes = [g for g in genes if g in panel]
    # Representative only: all directly affected genes remain explicit.
    selected_genes = panel_genes or genes
    impact = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
    gene = min(selected_genes, key=lambda g: (min(impact.get(t.get('IMPACT'), 4)
                                                for t in direct if t['SYMBOL'] == g), g))
    choices = [t for t in direct if t['SYMBOL'] == gene]
    preference = panel.get(gene, '')
    tx = min(choices, key=lambda t: transcript_rank(t, preference))
    warnings = []
    if preference and not any(matches_preference(t, preference) for t in choices):
        warnings.append('Preferred transcript ' + preference + ' unavailable for ' + gene + '; used matching-gene fallback')
    return dict(transcript=tx, genes=genes, ambiguous=len(genes) > 1, warnings=warnings)


def display_transcript(rows, panel=None):
    return select_transcript(rows, panel)['transcript']
