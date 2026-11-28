#!/usr/bin/env python3
"""Prepare additive Classy blood report artifacts from unmodified native scores."""
import argparse
import csv
from html import escape
import json
from pathlib import Path
import shutil
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'review'))
from classy_blood import classifiers, native_rows, class_label, finite, LABELS


def prepare(input, sample, plots, output, enable_alma3, modification_mode, image):
    raw = json.loads(Path(input).read_text())
    if raw.get('sample_label', raw.get('sample_id', sample)) != sample:
        raise ValueError('Classifier sample mismatch')
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    raw['oncoseq_classifier'] = dict(workflow='classy_blood', enable_alma3=enable_alma3,
                                    alma3_modification_mode=modification_mode, container=image,
                                    source_digest='sha256:cd1f5ea7a1ecf140f7ff5da013018f5f5f22c257bf8ad861899d2882d6dc3e91')
    (root / 'review_classification.json').write_text(json.dumps(raw, allow_nan=False)+'\n')
    summaries = classifiers(raw)
    tables = root / 'tables'
    tables.mkdir(exist_ok=True)
    content = ['# Methylation', '', 'Blood-focused classification by Classy. Native model scores, cutoffs and accepted-call decisions are retained.', '']
    if enable_alma3:
        content += ['Analysis performed using ALMA3. Research use only. ALMA3.1 is a native approximation; numerical parity with the reference implementation is not established.', '']
    inputs = root / 'native_report_inputs'
    inputs.mkdir(exist_ok=True)
    for plot in sorted(Path(plots).glob('*.svg')):
        if plot.name.endswith('_combined_blood.svg'):
            shutil.copyfile(plot, inputs / plot.name)
            content += ['## Blood classifiers', '', f'![]({plot.name})', '']
    groups = native_rows(raw)
    for summary in summaries:
        key = summary['name']
        entries = sorted(groups[key], key=lambda r: -(r.get('score') if finite(r.get('score')) else -1))[:5]
        if key == 'alma3':
            entries = []
            for task in sorted({r['classifier_run_id'] for r in groups[key]}):
                entries.extend(sorted([r for r in groups[key] if r['classifier_run_id'] == task],
                                      key=lambda r: -(r.get('score') if finite(r.get('score')) else -1))[:3])
        rows = []
        for row in entries:
            score, cutoff = row.get('score'), row.get('cutoff')
            accepted = row.get('accepted_call')
            decision = ('Unavailable' if summary['status'] != 'completed' else
                        'Accepted' if accepted is True else 'Not accepted' if accepted is False else
                        'Above cutoff' if finite(score) and finite(cutoff) and score >= cutoff else
                        'Below cutoff' if finite(score) and finite(cutoff) else 'Not reported')
            rows.append([row['classifier_run_id'], class_label(row), score if finite(score) and summary['status'] == 'completed' else 'NA', cutoff if finite(cutoff) else 'NA', decision, summary['status']])
        if not rows:
            rows = [[key, summary['reason'] or 'No result', 'NA', 'NA', 'Unavailable', summary['status']]]
        headers = ['Native task', 'Class / hierarchy', 'Score', 'Native cutoff', 'Native decision', 'Status']
        with (tables / (key+'.tsv')).open('w') as out:
            writer = csv.writer(out, delimiter='\t');writer.writerow(headers);writer.writerows(rows)
        content += ['## '+summary['display_name'], '', summary['reason'] or '', '',
                    '<table><thead><tr>'+''.join('<th>'+escape(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+escape(str(v))+'</td>' for v in r)+'</tr>' for r in rows)+'</tbody></table>', '']
    (inputs / (sample+'-Methylation.qmd')).write_text('\n'.join(content))
    metadata = dict(sample_id=sample, source='classy_blood', settings=raw['oncoseq_classifier'],
                    plots=[p.name for p in sorted(Path(plots).glob('*.svg'))],
                    models=[{k:r.get(k) for k in ('name','status','reason','model','model_version','model_version_source','target_coverage')} for r in summaries])
    (root / 'classification.metadata.json').write_text(json.dumps(metadata, allow_nan=False)+'\n')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('input','sample','plots','output','image'):p.add_argument('--'+name,required=True)
    p.add_argument('--enable-alma3',choices=['true','false'],default='true')
    p.add_argument('--modification-mode',choices=['5mc','5mc-plus-5hmc'],default='5mc')
    a=p.parse_args();prepare(a.input,a.sample,a.plots,a.output,a.enable_alma3=='true',a.modification_mode,a.image)
