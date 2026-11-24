"""Build versioned review data from native analysis artifacts."""
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import pysam

from adapters import classifiers, coverage, interval_table, nasvar_events
from annotate_small_variants import prepare as prepare_catalog, finish as finish_catalog
from annotate_structural_variants import annotate as annotate_sv
from report_prioritization import fusion_pairs, POLICY_VERSION
from resources import validate
from structural_consensus import consolidate
from transcript_selection import panel_provenance
from progress import progress


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def prepare(sample, inputs, references, panel, output, payload):
    refs, manifest = validate(references)
    entries = json.loads(Path(inputs).read_text())
    out, saved = Path(output), Path(payload)
    out.mkdir(parents=True, exist_ok=True)
    saved.mkdir(parents=True, exist_ok=True)
    sources, inventory = [], []
    for ordinal, entry in enumerate(entries):
        source = Path(entry['path'])
        if not source.is_file():
            raise ValueError('Missing review artifact: ' + str(source))
        target = saved / f'{ordinal:04d}_{source.name}'
        shutil.copyfile(source, target)
        item = dict(entry, path=target.name)
        inventory.append(item)
        if entry['kind'].startswith('small:'):
            caller = entry['kind'].split(':')[1]
            sources.append(dict(path=str(source.resolve()), sample=sample, caller=caller,
                                kind='germline' if caller == 'clair3' else 'tumor_only'))
    write_json(saved / 'inputs.json', inventory)
    write_json(saved / 'references.json', manifest)
    write_json(saved / 'panel.json', panel_provenance(panel))
    write_json(out / 'native_sources.json', sources)
    if sources:
        prepare_catalog(out / 'native_sources.json', refs['fasta'], refs['targets'], out)
        from allele_consensus import add
        add(out)
    else:
        write_json(out / 'summary.json', dict(candidate_alleles=0, status='not_run'))
    header=pysam.VariantHeader()
    header.info.add('CSQ','.', 'String','Consequence annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Feature|BIOTYPE|MANE_SELECT|MANE_PLUS_CLINICAL|CANONICAL|HGVSc|HGVSp|FLAGS')
    with pysam.VariantFile(out/'empty.vep.vcf','w',header=header):pass
    return refs


def complete(sample, catalogue, payload, vep, references, panel, known_fusions, distance, output, image_lock=None, code_digest=None):
    csv.field_size_limit(100_000_000)
    progress('Review completion started; validating reference files')
    refs, manifest = validate(references)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    cat, saved = Path(catalogue), Path(payload)
    inventory = json.loads((saved / 'inputs.json').read_text())
    sources = json.loads((cat / 'native_sources.json').read_text())
    if sources:
        finish_catalog(cat, vep, refs, sample, out, panel)
        small = json.loads((out / 'annotation' / (sample + '.small_variant_review.json')).read_text())
    else:
        small = dict(status='not_run', records=[], total_records=0)
    events, classification, cov, cnv, native, figures = [], [], [], [], {}, []
    import sqlite3
    from breakpoint_annotation import gene_index
    gene_path=out/'.genes.sqlite'
    gene_db=sqlite3.connect(gene_path)
    progress('Building structural gene index')
    gene_audit=gene_index(refs['gff'],refs['fai'],gene_db)
    gene_db.close()
    statuses = {}
    for item in inventory:
        kind, path = item['kind'], saved / item['path']
        if kind.startswith('small:'):
            statuses[kind] = 'completed'
        elif kind.startswith('sv:'):
            caller = kind.split(':')[1]
            folder = out / caller
            progress('Annotating structural source: ' + caller)
            annotate_sv(sample, caller, path, refs['gff'], refs['fai'], folder, (gene_path,gene_audit))
            for f in folder.glob('*.evidence.tsv.gz'):
                with gzip.open(f, 'rt') as stream:
                    events.extend(json.loads(r['event_json']) for r in csv.DictReader(stream, delimiter='\t'))
            statuses[kind] = 'completed'
        elif kind == 'classification':
            raw = json.loads(path.read_text())
            if raw.get('sample_label', raw.get('sample_id', sample)) != sample:
                raise ValueError('Classifier sample mismatch')
            classification = classifiers(raw)
            shutil.copyfile(path, out / 'classification.native.json')
            statuses[kind] = 'completed'
        elif kind == 'nasvar':
            native = json.loads(path.read_text())
            events.extend(nasvar_events(native, sample, (gene_path, gene_audit)))
            shutil.copyfile(path, out / 'nasvar.native.json')
            statuses[kind] = 'completed'
        elif kind == 'qc:alignment':
            shutil.copyfile(path, out / 'qc.alignment.json')
        elif kind in ('qc:cramino', 'qc:seqkit'):
            shutil.copyfile(path, out / ('qc.' + kind.split(':')[1] + '.tsv'))
        elif kind in ('nasvar:coverage', 'nasvar:maf'):
            shutil.copyfile(path, out / ('nasvar.' + kind.split(':')[1] + '.tsv'))
        elif kind.startswith('coverage:'):
            cov.append((kind, path))
        elif kind.startswith('cnv:'):
            cnv.extend(dict(r, source=kind.split(':')[1], coordinate_convention='BED 0-based start, end-exclusive' if kind=='cnv:qdnaseq' else 'Native VCF 1-based POS, inclusive END') for r in interval_table(path))
        elif kind.startswith('figure:'):
            folder=out/'figures'
            folder.mkdir(exist_ok=True)
            shutil.copyfile(path,folder/path.name)
            figures.append(dict(category=kind.split(':')[1], path='figures/'+path.name,
                                name=item.get('name', path.name)))
    gene_path.unlink()
    progress(f'Consolidating structural junctions: {len(events)} records')
    structural = consolidate(events, fusion_pairs(known_fusions), distance)
    progress(f'Structural consolidation finished: {len(structural)} events')
    for name in ('small:clair3', 'small:clairsto', 'small:deepsomatic',
                 'sv:sniffles', 'sv:severus', 'sv:stellerator', 'classification', 'nasvar'):
        statuses.setdefault(name, 'not_run')
    images=json.loads(Path(image_lock).read_text()) if image_lock else {}
    result = dict(schema_version=1, sample=sample, assembly='GRCh38', mode='review', runtime_version='1.0.0', code_sha256=code_digest, images_provenance=images,
                  policy_version=POLICY_VERSION, fusion_distance=distance,
                  statuses=statuses, small_variants=small, structural_events=structural,
                  classification=classification, coverage=coverage(cov), cnv=cnv,
                  nasvar=native, figures=figures, references=manifest,
                  panel=panel_provenance(panel), warnings=[])
    if manifest.get('pinning') == 'development':
        result['warnings'].append('Development references are unpinned; reference checksums must be pinned before clinical launch.')
    progress('Writing final review data and provenance')
    write_json(out / 'review.json', result)
    write_json(out / 'manifest.json', dict(schema_version=1, sample=sample, assembly='GRCh38',
                                         policy_version=POLICY_VERSION, references=manifest, images=images, code_sha256=code_digest,
                                         panel=result['panel'], statuses=statuses,
                                         inputs=[{k: v for k, v in i.items() if k != 'path'} for i in inventory]))
    with gzip.open(out / 'structural_events.tsv.gz', 'wt') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['event_id', 'selected', 'consensus', 'event_json'])
        for event in structural:
            writer.writerow([event['id'], event['selected'], event['consensus'], json.dumps(event)])
    with gzip.open(out / 'small_variants.tsv.gz', 'wt') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['allele', 'caller', 'kind', 'filter', 'population_status', 'annotation', 'source_record', 'native_records'])
        for path in sorted(out.glob('*/*.annotation.tsv.gz')):
            with gzip.open(path, 'rt') as source:
                reader = csv.reader(source, delimiter='\t')
                next(reader)
                writer.writerows(reader)
    with gzip.open(out / 'caller_evidence.tsv.gz', 'wt') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['allele', 'caller_evidence_json'])
        for row in small['records']:
            writer.writerow([json.dumps(row['coordinates']), json.dumps(row['evidence'])])
    progress('Review completion finished')
    return result


def import_results(results, sample, output):
    """Create a strict native-artifact inventory for a completed execution."""
    root, out = Path(results), Path(output)
    if not root.is_dir():
        raise ValueError('Results directory does not exist')
    inventory = []
    patterns = {
        'small:clair3': f'program_outputs/clair3/{sample}/merge_output.vcf.gz',
        'small:clairsto': f'program_outputs/bcftools/{sample}_snp_clairsto.vcf.gz',
        'small:deepsomatic': f'program_outputs/deepsomatic/{sample}.vcf.gz',
        'sv:sniffles': f'program_outputs/sniffles/{sample}_sv_sniffles.vcf.gz',
        'sv:severus': f'program_outputs/htslib/{sample}_sv_severus.vcf.gz',
        'sv:stellerator': f'program_outputs/htslib/{sample}_sv_stellerator.vcf.gz',
        'classification': f'reports/tumor_classifiers/{sample}_combined_classification.json',
        'nasvar': f'nasvar/{sample}.result.json',
        'coverage:primary': f'reports/mosdepth/{sample}_panel_primary.regions.bed.gz',
        'coverage:nofilter': f'program_outputs/mosdepth/{sample}_panel_nofilter.regions.bed.gz',
        'coverage:mapq60': f'program_outputs/mosdepth/{sample}_panel_mapq60.regions.bed.gz',
    }
    missing = []
    for kind, pattern in patterns.items():
        path = root / pattern
        if not path.is_file():
            missing.append(kind)
        else:
            inventory.append(dict(kind=kind, path=str(path.resolve()), name=path.name))
    if missing:
        raise ValueError('Missing native artifacts: ' + ', '.join(missing))
    for path in sorted((root / 'program_outputs' / 'figeno').glob(sample + '*')):
        if path.suffix.lower() in ('.png', '.svg'):
            category = 'fusions' if 'fusion' in path.name.lower() else 'target_loci' if 'target' in path.name.lower() else 'chromosomes'
            inventory.append(dict(kind='figure:' + category, path=str(path.resolve()), name=path.name))
    out.parent.mkdir(parents=True, exist_ok=True)
    write_json(out, inventory)
    return inventory
