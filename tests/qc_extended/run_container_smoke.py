#!/usr/bin/env python3
"""Execute the locked Linux amd64 QC tools on synthetic modification data."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

METRICS_IMAGE = 'community.wave.seqera.io/library/oncoseq_review_1_0_0@sha256:b143e350e1dae33ad26758fc6f9ccf2392cb8b49881de73960f42525602d5eb5'
MODKIT_IMAGE = 'ghcr.io/chusj-pigu/classy-slim@sha256:7ec5fc1e4607be13d30dc65a1b07f2844f8289edbec0ca67a70e9b88ef8beb4c'

parser = argparse.ArgumentParser()
parser.add_argument('--threads', type=int, default=4)
a = parser.parse_args()
if a.threads < 1:
    parser.error('threads must be positive')
repo = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='synthetic-qc-container-') as tmp:
    def run(image, script):
        subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64', '--network', 'none',
                        '--entrypoint', '/bin/bash', '-v', tmp + ':/data', '-v', str(repo / 'bin') + ':/workflow/bin:ro',
                        '-w', '/data', image, '-euc', ("export PATH=/opt/wave/.pixi/envs/default/bin:$PATH\n" + script if image == METRICS_IMAGE else script)], check=True)
    run(METRICS_IMAGE, '''python3 - <<'PY'
from pathlib import Path
from array import array
import pysam
Path('reference.fa').write_text('>chr1\\n'+'CG'*1000+'\\n')
pysam.faidx('reference.fa')
Path('enrichment.bed').write_text('chr1\\t0\\t1000\\n')
Path('targets.bed').write_text('chr1\\t200\\t600\\tSYNTHETIC_TARGET\\n')
with pysam.AlignmentFile('tagged.bam','wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000}]}) as out:
    for i in range(1000):
        r=pysam.AlignedSegment(); r.query_name=f'SYNTHETIC_READ_{i}'; r.query_sequence='CG'*250
        r.flag=0; r.reference_id=0; r.reference_start=0; r.mapping_quality=60; r.cigarstring='500M'
        r.set_tag('MM','C+m,'+','.join(['0']*250)+';'); r.set_tag('ML',array('B',[240]*250)); r.set_tag('MN',500)
        out.write(r)
pysam.index('tagged.bam')
PY
/workflow/bin/extended_qc.py reference --fasta reference.fa --cpgs reference.cpgs.gz --output reference.json
/workflow/bin/extended_qc.py alignment --bam tagged.bam --reference reference.json --enrichment enrichment.bed --targets targets.bed --threads 2 --output alignment.json
''')
    run(MODKIT_IMAGE, f'''modkit pileup tagged.bam pileup.bed.gz --ref reference.fa --cpg --modified-bases 5mC 5hmC --combine-strands --bgzf --threads {a.threads} --num-reads 10042 --filter-percentile 0.1 --log-filepath modkit.log''')
    run(METRICS_IMAGE, '''/workflow/bin/extended_qc.py methylation --alignment alignment.json --reference reference.json --cpgs reference.cpgs.gz --fasta reference.fa --bedmethyl pileup.bed.gz --enrichment enrichment.bed --targets targets.bed --output methylation.json
/workflow/bin/extended_qc.py assemble --sample SYNTHETIC_SMOKE --alignment alignment.json --methylation methylation.json --reference reference.json --enrichment enrichment.bed --targets targets.bed --output metrics
''')
    data = json.loads((Path(tmp) / 'metrics/metrics.json').read_text())
    primary = next(r for r in data['metrics'] if r['metric_id'] == 'alignment.primary_records')
    sites = next(r for r in data['metrics'] if r['metric_id'] == 'methylation.sites' and r['scope'] == 'genome_primary' and r['modification'] == 'combined')
    assert primary['value'] == 1000 and sites['value'] > 0
    assert data['components']['methylation']['status'] == 'complete'
print('Actual locked-container modkit and QC smoke: PASS')
