"""Recorded run provenance; configured values never masquerade as observations."""
import csv
import json
import re
from pathlib import Path

PARAMETERS = ('run_id', 'genome', 'basecall_model', 'm_bases', 'skip_basecalling',
              'skip_mapping', 'basecall_shards', 'device', 'bed', 'padding',
              'qc_targets_bed', 'time_series', 'cfdna', 'clinical', 'classy_blood',
              'classy_blood_enable_alma3', 'classy_blood_alma3_modification_mode')
WORKFLOW_FIELDS = ('pipeline', 'pipeline_version', 'nextflow_version', 'commit_id', 'revision')


def recorded_context(raw):
    return dict(workflow={k: v for k, v in raw.get('workflow', {}).items() if k in WORKFLOW_FIELDS},
                parameters={k: v for k, v in raw.get('parameters', {}).items() if k in PARAMETERS},
                software_version_texts=list(raw.get('software_version_texts', [])))


def version_rows(text, source):
    """Read the flat process/tool scalar format emitted by versions.yml.

    Keep scalars as strings so versions such as 1.20 do not become 1.2.
    No general YAML deserialization, tags, or object construction is used.
    """
    rows, process = [], None
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith(('#', '```')):
            continue
        if not line[0].isspace() and line.rstrip().endswith(':'):
            process = line.strip()[:-1].strip('"\'').rsplit(':', 1)[-1]
        elif process and line[0].isspace() and ':' in line:
            tool, value = line.strip().split(':', 1)
            value = value.strip()
            if value.startswith('"') and value.endswith('"'):
                try:
                    value = json.loads(value)
                except ValueError:
                    value = value[1:-1]
            elif value.startswith("'") and value.endswith("'"):
                value = value[1:-1].replace("''", "'")
            if value and value.lower() not in ('null', 'none', '~'):
                rows.append(dict(process=process, tool=tool.strip('"\''), version=value, source=source))
    return rows


def completed_run(root, sample):
    root = Path(root)
    context = dict(workflow={}, parameters={}, versions=[], containers=[])
    params = sorted((root / 'pipeline_info').glob('params_*.json'))
    if params:
        values = json.loads(params[-1].read_text())
        context['parameters'] = {k: values[k] for k in PARAMETERS if k in values}
    versions = root / 'reports' / (sample + '_report_output') / (sample + '-Versions-versions.qmd')
    if versions.is_file():
        context['versions'] = version_rows(versions.read_text(), 'Native workflow version report')
    observed = {}
    for path in sorted((root / 'pipeline_info').glob('execution_trace_*.txt')):
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if row.get('status') not in ('COMPLETED', 'CACHED') or not row.get('container') or row['container'] == '-':
                    continue
                process = row['name'].split(' (', 1)[0].rsplit(':', 1)[-1]
                key = (process, row['container'])
                observed[key] = dict(process=process, image=row['container'], source='Execution trace')
    context['containers'] = sorted(observed.values(), key=lambda r: (r['process'], r['image']))
    for row in context['versions']:
        if row['process'] == 'Workflow':
            context['workflow'][row['tool']] = row['version']
    return context


def enrich(data, results_root=None):
    recorded = data.get('runtime_provenance', {})
    context = dict(workflow=dict(recorded.get('workflow', {})),
                   parameters={k: v for k, v in recorded.get('parameters', {}).items() if k in PARAMETERS},
                   versions=[], containers=[])
    for text in recorded.get('software_version_texts', []):
        context['versions'].extend(version_rows(text, 'Captured versions.yml'))
    if results_root and data.get('sample'):
        completed = completed_run(results_root, data['sample'])
        context['workflow'].update(completed['workflow'])
        context['parameters'].update(completed['parameters'])
        context['versions'].extend(completed['versions'])
        context['containers'] = completed['containers']
    qc = data.get('extended_qc', {}).get('provenance', {})
    for tool, version in qc.get('versions', {}).items():
        if tool != 'depth_backend':
            context['versions'].append(dict(process='Extended QC alignment', tool=tool, version=str(version), source='Recorded QC metrics'))
    if qc.get('modkit_version') not in (None, 'not_run'):
        context['versions'].append(dict(process='Extended QC methylation', tool='modkit', version=str(qc['modkit_version']), source='Recorded QC metrics'))
    unique = {(r['process'], r['tool'], str(r['version']), r['source']): r for r in context['versions']}
    context['versions'] = sorted(unique.values(), key=lambda r: (r['tool'].lower(), r['process'], str(r['version'])))
    models = []
    for model in qc.get('basecalling_models', []):
        models.append(dict(kind='Basecalling / modified-base model', specification=model, source='BAM read-group metadata'))
    params = context['parameters']
    if params.get('basecall_model'):
        models.append(dict(kind='Basecalling model', specification=params['basecall_model'], source='Configured run parameters'))
    if params.get('m_bases'):
        models.append(dict(kind='Modified-base model request', specification=params['m_bases'], source='Configured run parameters'))
        if params.get('basecall_model') and '@' in params['basecall_model'] and '@' in params['m_bases']:
            models.append(dict(kind='Modified-base model specification', specification=params['basecall_model'] + '_' + params['m_bases'], source='Derived from configured model requests'))
    for classifier in data.get('classification', []):
        if classifier.get('model_version'):
            models.append(dict(kind=classifier.get('display_name', classifier['name']) + ' classifier',
                               specification=str(classifier['model_version']),
                               source=classifier.get('model_version_source') or 'Native classifier report'))
    context['models'] = models
    classifier = data.get('classifier_provenance', {})
    if classifier.get('workflow') == 'classy_blood' and classifier.get('container'):
        context['containers'].append(dict(process='CLASSY_BLOOD', image=classifier['container'],
                                          source='Captured classifier invocation'))
    lock = data.get('images_provenance', {})
    context['configured_containers'] = [dict(name=name, image=entry['image'], processes=entry.get('processes', [])) for name, entry in lock.get('images', {}).items()]
    context['qc_containers'] = [dict(process=process, image=image, source='Recorded QC metric provenance') for process, image in qc.get('images', {}).items()]
    context['identity'] = dict(assembly=data.get('assembly'), review_runtime=data.get('runtime_version'),
                               review_code_sha256=data.get('code_sha256'), policy_version=data.get('policy_version'),
                               container_platform=lock.get('platform'))
    context['checksums'] = [dict(asset=key.removesuffix('_sha256').replace('_', ' ').capitalize(), sha256=qc[key], source='Recorded QC metrics')
                            for key in ('reference_sha256', 'targets_sha256', 'enrichment_sha256') if qc.get(key)]
    if classifier.get('source_digest'):
        context['checksums'].append(dict(asset='Classy blood source image', sha256=classifier['source_digest'].removeprefix('sha256:'),
                                         source='Captured classifier invocation'))
    for key, entry in data.get('references', {}).get('resources', {}).items():
        if isinstance(entry, dict) and entry.get('sha256'):
            context['checksums'].append(dict(asset=key, sha256=entry['sha256'], source='Review reference manifest'))
    data['run_provenance'] = context


def attach_report(input, output, versions):
    source = Path(input).read_text()
    match = re.search(r'(<script type="application/json" id="review-data">)(.*?)(</script>)', source, re.S)
    if not match:
        raise ValueError('Review report data block not found')
    data = json.loads(match.group(2))
    previous = data.get('run_provenance', {})
    context = recorded_context(data.get('runtime_provenance', {}))
    context['software_version_texts'] = json.loads(Path(versions).read_text())
    data['runtime_provenance'] = context
    enrich(data)
    # Retain execution-trace observations when attaching to a local rerender.
    if previous.get('containers'):
        data['run_provenance']['containers'] = previous['containers']
    payload = json.dumps(data, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c').replace('&', '\\u0026')
    Path(output).write_text(source[:match.start(2)] + payload + source[match.end(2):])
