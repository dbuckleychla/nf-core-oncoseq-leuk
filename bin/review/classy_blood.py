"""Adapt Classy blood reports without changing scores or native decisions."""
import math
from collections import defaultdict

LABELS = {'marlin': 'MARLIN', 'lamprey': 'Lamprey', 'alma3': 'ALMA3.1 (native approximation)', 'optic': 'OPTiC'}


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def family(run):
    name = str(run).lower()
    return next((key for key in LABELS if name == key or name.startswith(key + ':') or name.startswith(key + '_') or name.startswith(key + '.')), None)


def native_rows(data):
    grouped = defaultdict(list)
    for row in data.get('combined', {}).get('classification_scores', []):
        key = family(row.get('classifier_run_id'))
        if key:
            grouped[key].append(row)
    return grouped


def class_label(row):
    levels = [str(p['name']) for p in row.get('native_path', []) if isinstance(p, dict) and p.get('name')]
    return ' / '.join(levels) or str(row.get('class_name') or row.get('name') or row.get('label') or 'Unreported class')


def classifiers(data):
    groups = native_rows(data)
    settings = data.get('oncoseq_classifier', {})
    result = []
    for key, name in LABELS.items():
        model = data.get(key) or {}
        rows = sorted(groups[key], key=lambda r: -(r.get('score') if finite(r.get('score')) else -1))
        coverage = model.get('target_coverage') or {}
        observed = coverage.get('observed_target_count', model.get('processed_observations'))
        status, reason = 'completed', None
        if key == 'alma3' and settings.get('enable_alma3') is False:
            status, reason = 'disabled', 'ALMA3 native approximation disabled by configuration'
        elif not rows:
            status, reason = 'unavailable', 'No native classification scores emitted'
        elif observed == 0:
            status, reason = 'unavailable', 'No matched CpG observations; emitted probabilities are unsupported'
        accepted = [r for r in rows if r.get('accepted_call') is True]
        top = (accepted[0] if accepted else rows[0]) if rows else {}
        if key == 'alma3':
            path = [p for p in model.get('result', {}).get('path', []) if p.get('status') == 'resolved']
            if path:
                native = path[-1]
                matches = [r for r in rows if r.get('classifier_run_id') == 'alma3:' + native.get('level', '')
                           and r.get('native_label') == native.get('classification')]
                if matches:
                    top = matches[0]
        cutoff = top.get('cutoff', model.get('cutoff'))
        cutoff = cutoff if finite(cutoff) else None
        decision = top.get('accepted_call') if isinstance(top.get('accepted_call'), bool) else None
        if status != 'completed':
            decision = None
        scores = [dict(label=class_label(r), score=r['score'], accepted_call=r.get('accepted_call')) for r in rows if finite(r.get('score'))]
        result.append(dict(name=key, display_name=name, status=status, reason=reason,
                           label=((model.get('result', {}).get('result_summary') if key == 'alma3' else None) or class_label(top)) if top and status == 'completed' else None,
                           score=top.get('score') if status == 'completed' and finite(top.get('score')) else None,
                           cutoff=cutoff, accepted_call=decision, thresholds={'cutoff': cutoff} if cutoff is not None else {},
                           scores=scores, interpretation='Native classification score', model=model.get('model_label'),
                           model_version=model.get('model_version'), model_version_source=model.get('model_version_source'),
                           target_coverage=coverage, native_report=model,
                           approximation=key == 'alma3'))
    return result
