"""Render a shareable layout preview from synthetic records and illustrations."""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin/review'))

from fixtures import make
from engine import prepare, complete
from render import render
from structural_consensus import consolidate


def preview(output):
    project = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory() as tmp:
        root = make(Path(tmp) / 'references')
        entries = json.loads((root / 'inputs.json').read_text())
        for entry in entries:
            entry['path'] = str(root / entry['path'])
        (root / 'inputs.json').write_text(json.dumps(entries))
        prepare('SYNTHETIC', root / 'inputs.json', root / 'references.json', root / 'panel.tsv',
                Path(tmp) / 'catalogue', Path(tmp) / 'payload')
        data = complete('SYNTHETIC', Path(tmp) / 'catalogue', Path(tmp) / 'payload', root / 'fixture.vep.vcf.gz',
                        root / 'references.json', root / 'panel.tsv', project / 'assets/fusion-list.txt',
                        100, Path(tmp) / 'review', image_lock=project / 'images.lock.json')
        rows = []
        for caller, delta, reads in [('sniffles', 0, 65), ('severus', 5, 63)]:
            rows.append(dict(id='illustrative:' + caller, caller=caller, sample='SYNTHETIC', chrom='chr1', pos=800 + delta,
                             partner=dict(chrom='chr2', pos=900 + delta, local_side='right', remote_side='left'),
                             filter=['PASS'], evidence=dict(info=dict(SVTYPE='BND', SUPPORT=reads), format={}),
                             genes1=[dict(name='RUNX1')], genes2=[dict(name='RUNX1T1')],
                             gene_pair=['RUNX1', 'RUNX1T1'], source_records=[], coordinate_convention='Illustrative coordinates'))
        data['structural_events'] = consolidate(rows, {('RUNX1', 'RUNX1T1')})
        data['mode'] = 'preview'
        data['warnings'] = ['Synthetic records and an illustrative fusion figure are used to demonstrate the layout.']
        figures = Path(tmp) / 'figures'
        figures.mkdir()
        (figures / 'illustrative.svg').write_text('''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="440" viewBox="0 0 1100 440">
<rect width="1100" height="440" fill="white"/><text x="35" y="38" font-family="sans-serif" font-size="22" fill="#172c3d">Illustrative junction: RUNX1 :: RUNX1T1</text>
<text x="35" y="76" font-family="sans-serif" font-size="18" fill="#526774">Synthetic coverage and split-read tracks</text>
<g stroke="#d7e0e6"><path d="M120 110H510M590 110H1010M120 220H510M590 220H1010M120 350H510M590 350H1010"/></g>
<g font-family="sans-serif" font-size="18" fill="#172c3d"><text x="15" y="150">Coverage</text><text x="15" y="270">Split reads</text><text x="15" y="395">Genes</text><text x="290" y="420">RUNX1</text><text x="760" y="420">RUNX1T1</text></g>
<g fill="#748a96"><path d="M130 210V155H160V150H190V140H220V130H250V135H280V140H310V130H340V125H370V145H400V130H430V135H460V145H490V210Z"/><path d="M610 210V145H640V130H670V140H700V125H730V135H760V130H790V120H820V135H850V145H880V140H910V135H940V145H970V150H1000V210Z"/></g>
<g stroke="#086f83" stroke-width="6"><path d="M170 260H410M160 280H420M180 300H400M700 260H940M710 280H950M700 300H930"/></g>
<g stroke="#955715" stroke-width="2" fill="none"><path d="M410 260Q555 220 700 260M420 280Q550 240 710 280M400 300Q550 265 700 300"/></g>
<g fill="#086f83"><rect x="200" y="375" width="55" height="20"/><rect x="300" y="375" width="40" height="20"/><rect x="405" y="375" width="60" height="20"/></g>
<g fill="#b16832"><rect x="630" y="375" width="45" height="20"/><rect x="745" y="375" width="60" height="20"/><rect x="905" y="375" width="60" height="20"/></g>
</svg>''')
        manifest = [dict(category='fusions', generated=True, event_id=data['structural_events'][0]['id'],
                         path='illustrative.svg', name='Illustrative RUNX1::RUNX1T1 junction')]
        (figures / 'figures.json').write_text(json.dumps(manifest))
        data_path = Path(tmp) / 'review/review.json'
        data_path.write_text(json.dumps(data))
        render(data_path, output, project / 'assets/review/report.html', figures)


if __name__ == '__main__':
    preview(sys.argv[1])
