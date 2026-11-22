"""Small synthetic references and native records for review regression tests."""
import hashlib
import json
from array import array
from pathlib import Path

import pysam


def make(directory):
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / 'genome.fa'
    fasta.write_text('>chr1\n' + 'A' * 1000 + '\n>chr2\n' + 'A' * 1000 + '\n')
    pysam.faidx(str(fasta))
    (root / 'targets.bed').write_text('chr1\t0\t1000\tTP53\n')
    (root / 'genes.gff').write_text('##gff-version 3\nchr1\ttest\tgene\t100\t300\t.\t+\t.\tID=gene-TP53;Name=TP53\nchr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=NM_TEST.1;Parent=gene-TP53\nchr1\ttest\texon\t100\t300\t.\t+\t.\tParent=NM_TEST.1\n')
    (root / 'panel.tsv').write_text('gene\tpreferred_transcript\nTP53\tNM_TEST.1\n')

    def header(sample=False):
        h = pysam.VariantHeader()
        for chrom in ('chr1', 'chr2'):
            h.contigs.add(chrom, length=1000)
        if sample:
            h.add_sample('SYNTHETIC')
            for name, num, typ in [('GT', 1, 'String'), ('DP', 1, 'Integer'), ('AD', 'R', 'Integer'), ('AF', 'A', 'Float')]:
                h.formats.add(name, num, typ, 'Synthetic fixture')
        return h

    def write(name, h, rows):
        dest = root / name
        with pysam.VariantFile(dest, 'wz', header=h) as out:
            for pos, alt, identifier, info, passing in rows:
                rec = out.new_record(contig='chr1', start=pos - 1, alleles=('A', alt), id=identifier)
                for k, v in info.items():
                    rec.info[k] = v
                if passing:
                    rec.filter.add('PASS')
                if rec.samples:
                    rec.samples['SYNTHETIC']['GT'] = (0, 1)
                    rec.samples['SYNTHETIC']['DP'] = 100
                    rec.samples['SYNTHETIC']['AD'] = (70, 30)
                    rec.samples['SYNTHETIC']['AF'] = (0.3,)
                out.write(rec)
        pysam.tabix_index(str(dest), preset='vcf', force=True)
        return dest

    calls = [(150, 'T', 'call1', {}, True), (160, 'AT', 'call2', {}, True), (170, 'C', 'call3', {}, True)]
    for caller in ('clair3', 'clairsto', 'deepsomatic'):
        write(caller + '.vcf.gz', header(True), calls)
    h = header()
    h.info.add('AC', 'A', 'Integer', 'Synthetic allele count')
    h.info.add('AN', 1, 'Integer', 'Synthetic allele number')
    for name in ('gnomad_genome', 'gnomad_exome'):
        write(name + '.vcf.gz', h, [(150, 'T', '.', dict(AC=(1,), AN=1000), False),
                                  (160, 'AT', '.', dict(AC=(1,), AN=1000), False),
                                  (170, 'C', '.', dict(AC=(30,), AN=1000), False)])
    h = header()
    h.info.add('CLNSIG', '.', 'String', 'Synthetic significance')
    h.info.add('CLNREVSTAT', '.', 'String', 'Synthetic review status')
    write('clinvar.vcf.gz', h, [(150, 'T', '1', dict(CLNSIG=('Pathogenic',), CLNREVSTAT=('reviewed_by_expert_panel',)), False)])
    h = header()
    h.info.add('clinvar_cmn_path', 1, 'Integer', 'Explicit rescue')
    write('rescue.vcf.gz', h, [])
    h = header()
    for name, typ in [('GENE', 'String'), ('TRANSCRIPT', 'String'), ('TIER', 'String'), ('TARGETED_SCREEN_SAMPLE_COUNT', 'Integer')]:
        h.info.add(name, 1, typ, 'Synthetic COSMIC field')
    write('cosmic.vcf.gz', h, [(150, 'T', 'COSV1', dict(GENE='TP53', TRANSCRIPT='TX1', TIER='1', TARGETED_SCREEN_SAMPLE_COUNT=7), False),
                             (150, 'T', 'COSV1', dict(GENE='TP53', TRANSCRIPT='TX1', TIER='1', TARGETED_SCREEN_SAMPLE_COUNT=7), False),
                             (150, 'T', 'COSV1', dict(GENE='TP53', TRANSCRIPT='TX2', TIER='1', TARGETED_SCREEN_SAMPLE_COUNT=9), False)])
    h = header(True)
    for name, num, typ in [('SVTYPE', 1, 'String'), ('SUPPORT', 1, 'Integer'), ('MATEID', '.', 'String')]:
        h.info.add(name, num, typ, 'Synthetic SV field')
    for caller in ('sniffles','severus','stellerator'):
        write(caller+'.sv.vcf.gz', h, [(800, 'A[chr2:900[', 'sv1', dict(SVTYPE='BND', SUPPORT=5), True)])
    fields = 'Allele|Consequence|IMPACT|SYMBOL|Feature|BIOTYPE|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|HGVSc|HGVSp|FLAGS'
    h = header()
    h.info.add('CSQ', '.', 'String', 'Consequence annotations. Format: ' + fields)
    rows = []
    for pos, alt, consequence, impact in [(150, 'T', 'missense_variant', 'MODERATE'), (160, 'AT', 'frameshift_variant', 'HIGH'), (170, 'C', 'missense_variant', 'MODERATE')]:
        csq = '|'.join([alt, consequence, impact, 'TP53', 'NM_TEST.1', 'protein_coding', 'YES', 'NM_TEST.1', '', f'NM_TEST.1:c.{pos}A>T', 'NP_TEST:p.X', ''])
        rows.append((pos, alt, '.', dict(CSQ=(csq,)), False))
    write('fixture.vep.vcf.gz', h, rows)
    cache = root / 'cache'
    cache.mkdir(exist_ok=True)
    (cache / 'fixture.vep.vcf.gz').write_bytes((root / 'fixture.vep.vcf.gz').read_bytes())
    (root / 'classification.json').write_text(json.dumps(dict(sample_label='SYNTHETIC', marlin=dict(cutoff=0.8, inference=dict(annotated_probabilities=[dict(class_name='Synthetic model class', probability=0.9)])))))
    (root / 'nasvar.json').write_text(json.dumps(dict(version='synthetic', fusions=dict(fusions=[]), qc=dict(mean_depth=100))))
    for label, depth in [('primary', 80), ('mapq60', 60), ('nofilter', 100)]:
        (root / (label + '.bed')).write_text(f'chr1\t0\t1000\tTP53\t{depth}\n')
    bam = root / 'sample.bam'
    bh = dict(HD=dict(VN='1.6', SO='coordinate'), SQ=[dict(SN='chr1', LN=1000), dict(SN='chr2', LN=1000)])
    with pysam.AlignmentFile(bam, 'wb', header=bh) as out:
        for i, qs in enumerate((5, 15)):
            r = pysam.AlignedSegment(out.header)
            r.query_name = 'read' + str(i)
            r.query_sequence = 'AC' * 25
            r.query_qualities = array('B', [30] * 50)
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 750 + i
            r.mapping_quality = 60
            r.cigarstring = '50M'
            r.set_tag('qs', qs)
            r.set_tag('MM', 'C+m,0;')
            r.set_tag('ML', array('B', [128]))
            r.set_tag('MN', 50)
            out.write(r)
    pysam.index(str(bam))
    resources = {}
    mapping = dict(fasta='genome.fa', fai='genome.fa.fai', targets='targets.bed', gff='genes.gff',
                   gnomad_genome='gnomad_genome.vcf.gz', gnomad_exome='gnomad_exome.vcf.gz',
                   clinvar='clinvar.vcf.gz', clinvar_common_pathogenic='rescue.vcf.gz', cosmic_targeted='cosmic.vcf.gz')
    for name, value in mapping.items():
        resources[name] = dict(uri=value, release='synthetic', sha256=hashlib.sha256((root / value).read_bytes()).hexdigest())
        if value.endswith('.vcf.gz'):
            resources[name]['index'] = dict(uri=value + '.tbi', sha256=hashlib.sha256((root / (value + '.tbi')).read_bytes()).hexdigest())
    (root / 'references.json').write_text(json.dumps(dict(schema_version=1, assembly='GRCh38', resources=resources), indent=2) + '\n')
    inventory = []
    for caller in ('clair3', 'clairsto', 'deepsomatic'):
        inventory.append(dict(kind='small:' + caller, path=caller + '.vcf.gz'))
    for caller in ('sniffles', 'severus', 'stellerator'):
        inventory.append(dict(kind='sv:' + caller, path=caller+'.sv.vcf.gz'))
    (root / 'nasvar.coverage.tsv').write_text('chromosome\tbin_start\tbin_end\tn_reads\tn_reads_gc_adj\nchr1\t0\t1000000\t10\t12\n')
    (root / 'nasvar.maf').write_text('chrom\tpos\tref\talt\nchr1\t100\t15\t5\n')
    (root / 'cramino.txt').write_text('Number of reads\t1000\nN50\t15000\n')
    inventory += [dict(kind='nasvar:coverage', path='nasvar.coverage.tsv'), dict(kind='nasvar:maf', path='nasvar.maf'), dict(kind='qc:cramino', path='cramino.txt')]
    inventory += [dict(kind='classification', path='classification.json'), dict(kind='nasvar', path='nasvar.json')]
    inventory += [dict(kind='coverage:' + label, path=label + '.bed') for label in ('primary', 'mapq60', 'nofilter')]
    (root / 'inputs.json').write_text(json.dumps(inventory, indent=2) + '\n')
    return root


if __name__ == '__main__':
    import sys
    make(sys.argv[1])
