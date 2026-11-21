"""Validate immutable staged review references and prepare annotation VCFs."""
import hashlib
import json
import re
from pathlib import Path

import pysam
import pysam.bcftools

REQUIRED = ('fasta', 'fai', 'targets', 'gff', 'gnomad_genome', 'gnomad_exome',
            'clinvar', 'clinvar_common_pathogenic')


def checksum(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate(path):
    manifest = json.loads(Path(path).read_text())
    if manifest.get('schema_version') != 1 or manifest.get('assembly') != 'GRCh38':
        raise ValueError('Review references require schema_version 1 and assembly GRCh38')
    entries = manifest.get('resources', {})
    if set(REQUIRED) - entries.keys():
        raise ValueError('Missing review resources: ' + ', '.join(sorted(set(REQUIRED) - entries.keys())))
    result = {}
    for name, entry in entries.items():
        value = Path(entry.get('path') or entry.get('uri', ''))
        if not value.is_absolute():value=Path(path).parent/value
        if not value.is_file():
            raise ValueError('Review reference must be staged locally: ' + name)
        expected = entry.get('sha256', '')
        if not re.fullmatch('[0-9a-f]{64}', expected) or checksum(value) != expected:
            raise ValueError('Review reference checksum mismatch: ' + name)
        result[name] = str(value.resolve())
        if name.startswith(('gnomad', 'clinvar', 'cosmic')):
            index = entry.get('index', {})
            index_path = Path(index.get('path') or index.get('uri', ''))
            if not index_path.is_absolute():index_path=Path(path).parent/index_path
            if not index_path.is_file() or checksum(index_path) != index.get('sha256'):
                raise ValueError('Missing or invalid review VCF index: ' + name)
            with pysam.VariantFile(str(value), index_filename=str(index_path)) as handle:
                first = next(iter(handle.header.contigs), None)
                if first is not None:
                    next(handle.fetch(first), None)
    if Path(result['fai']).resolve() != Path(result['fasta'] + '.fai').resolve():
        raise ValueError('FASTA index must be adjacent to FASTA and named <fasta>.fai')
    return result, manifest


def prepare_vcf(source, fasta, output):
    """Rename only recognised primary contigs, then normalise, sort and index."""
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    lengths = {line.split('\t')[0] for line in Path(str(fasta) + '.fai').read_text().splitlines()}
    with pysam.VariantFile(source) as handle:
        mapping = {}
        for contig in handle.header.contigs:
            candidate = 'chr' + contig if contig in {str(i) for i in range(1, 23)} | {'X', 'Y'} else contig
            if contig in ('M', 'MT') and 'chrM' in lengths:
                candidate = 'chrM'
            if candidate not in lengths:
                raise ValueError('Annotation contig absent from FASTA: ' + contig)
            mapping[contig] = candidate
    names = out.with_name(out.name + '.contigs.tsv')
    names.write_text(''.join(f'{a}\t{b}\n' for a, b in mapping.items()))
    renamed = out.with_name(out.name + '.renamed.vcf.gz')
    normal = out.with_name(out.name + '.normal.vcf.gz')
    pysam.bcftools.annotate('--rename-chrs', str(names), '-Oz', '-o', str(renamed), source, catch_stdout=False)
    pysam.bcftools.norm('-f', fasta, '-c', 'e', '-m', '-any', '-Oz', '-o', str(normal), str(renamed), catch_stdout=False)
    pysam.bcftools.sort('-Oz', '-o', str(out), str(normal), catch_stdout=False)
    pysam.tabix_index(str(out), preset='vcf', force=True)
    renamed.unlink()
    normal.unlink()
    audit = dict(schema_version=1, source=str(source), source_sha256=checksum(source),
                 fasta_sha256=checksum(fasta), contig_mapping=mapping,
                 output=str(out), sha256=checksum(out), index_sha256=checksum(str(out) + '.tbi'))
    out.with_name(out.name + '.preparation.json').write_text(json.dumps(audit, indent=2) + '\n')
    return audit


def make_manifest(resources, output):
    """Record local prepared files with their deployment URIs and checksums."""
    entries = json.loads(Path(resources).read_text())
    if 'resources' in entries:
        entries = entries['resources']
    if set(REQUIRED) - entries.keys():
        raise ValueError('Resource mapping is missing required references')
    result = {}
    for name, value in entries.items():
        entry = dict(value) if isinstance(value, dict) else dict(path=value)
        path = Path(entry['path'])
        if not path.is_file():
            raise ValueError('Prepared resource does not exist: ' + str(path))
        uri = entry.get('uri', str(path.resolve()))
        result[name] = dict(uri=uri, release=entry.get('release', 'operator-supplied'), sha256=checksum(path))
        if name.startswith(('gnomad', 'clinvar', 'cosmic')):
            indexes = [Path(str(path) + suffix) for suffix in ('.tbi', '.csi')]
            index = next((p for p in indexes if p.is_file()), None)
            if index is None:
                raise ValueError('Prepared VCF has no adjacent index: ' + str(path))
            suffix = str(index)[len(str(path)):]
            result[name]['index'] = dict(uri=uri + suffix, sha256=checksum(index))
    manifest = dict(schema_version=1, assembly='GRCh38', resources=result)
    Path(output).write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest
