"""Prepare a reusable, disk-backed COSMIC targeted-screen evidence reference."""
import csv
import gzip
import json
from pathlib import Path
import re
import sqlite3

BLOOD_SITE = 'haematopoietic_and_lymphoid_tissue'
PHENOTYPE_FIELDS = ('COSMIC_PHENOTYPE_ID', 'PRIMARY_SITE', 'SITE_SUBTYPE_1',
                    'SITE_SUBTYPE_2', 'SITE_SUBTYPE_3', 'PRIMARY_HISTOLOGY',
                    'HISTOLOGY_SUBTYPE_1', 'HISTOLOGY_SUBTYPE_2', 'HISTOLOGY_SUBTYPE_3')
EVIDENCE_FIELDS = ('GENOMIC_MUTATION_ID', 'COSMIC_SAMPLE_ID', 'COSMIC_PHENOTYPE_ID',
                   'PUBMED_PMID', 'MUTATION_SOMATIC_STATUS', 'GENE_SYMBOL',
                   'TRANSCRIPT_ACCESSION', 'HGVSC', 'HGVSP')


def rows(path, required):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if set(required) - set(reader.fieldnames or []):
            raise ValueError('Missing COSMIC TSV columns: ' + ', '.join(sorted(set(required) - set(reader.fieldnames or []))))
        for row in reader:
            if None in row or any(row.get(k) is None for k in required):
                raise ValueError('Malformed COSMIC TSV row')
            yield {k: v.strip() for k, v in row.items()}


def prepare_evidence(mutations, classification, vcf, release, output):
    """Only positive screens join to COSV IDs; never retain SAMPLE_NAME."""
    from resources import checksum
    from progress import progress
    output = Path(output)
    if output.exists():
        raise ValueError('COSMIC evidence output already exists')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.building')
    if temporary.exists():
        raise ValueError('COSMIC evidence build already exists')
    db = sqlite3.connect(temporary)
    try:
        db.executescript('CREATE TABLE metadata (json TEXT); '
                         'CREATE TABLE phenotypes (id TEXT PRIMARY KEY, json TEXT); '
                         'CREATE TABLE evidence (' + ','.join(f'{k} TEXT NOT NULL' for k in EVIDENCE_FIELDS) + ', PRIMARY KEY (' + ','.join(EVIDENCE_FIELDS) + ')) WITHOUT ROWID;')
        for row in rows(classification, PHENOTYPE_FIELDS):
            identifier = row['COSMIC_PHENOTYPE_ID']
            encoded = json.dumps({k: row[k] for k in PHENOTYPE_FIELDS}, sort_keys=True)
            previous = db.execute('SELECT json FROM phenotypes WHERE id=?', (identifier,)).fetchone()
            if previous and previous[0] != encoded:
                raise ValueError('Conflicting COSMIC phenotype classification')
            db.execute('INSERT OR IGNORE INTO phenotypes VALUES (?,?)', (identifier, encoded))
        counts = dict(input_rows=0, positive_rows=0, negative_rows=0, positive_without_cosv=0)
        batch = []
        for row in rows(mutations, (*EVIDENCE_FIELDS, 'POSITIVE_SCREEN')):
            counts['input_rows'] += 1
            if counts['input_rows'] % 1000000 == 0:
                progress(f"COSMIC preparation: {counts['input_rows']} input rows")
            flag = row['POSITIVE_SCREEN'].lower()
            if flag not in ('y', 'n'):
                raise ValueError('Unexpected COSMIC POSITIVE_SCREEN value')
            if flag == 'n':
                counts['negative_rows'] += 1
                continue
            counts['positive_rows'] += 1
            if not row['GENOMIC_MUTATION_ID'] or row['GENOMIC_MUTATION_ID'] in ('NS', '.'):
                counts['positive_without_cosv'] += 1
                continue
            if not re.fullmatch(r'COSV\d+', row['GENOMIC_MUTATION_ID']):
                raise ValueError('Unexpected COSMIC genomic mutation identifier')
            batch.append(tuple(row[k] for k in EVIDENCE_FIELDS))
            if len(batch) == 10000:
                db.executemany('INSERT OR IGNORE INTO evidence VALUES (' + ','.join('?' for _ in EVIDENCE_FIELDS) + ')', batch)
                db.commit()
                batch.clear()
        db.executemany('INSERT OR IGNORE INTO evidence VALUES (' + ','.join('?' for _ in EVIDENCE_FIELDS) + ')', batch)
        counts['distinct_evidence_rows'] = db.execute('SELECT COUNT(*) FROM evidence').fetchone()[0]
        counts['unclassified_evidence_rows'] = db.execute('SELECT COUNT(*) FROM evidence e LEFT JOIN phenotypes p ON e.COSMIC_PHENOTYPE_ID=p.id WHERE p.id IS NULL').fetchone()[0]
        if not counts['distinct_evidence_rows']:
            raise ValueError('COSMIC targeted screens contain no positive COSV evidence')
        metadata = dict(schema_version=1, assembly='GRCh38', release=release, counts=counts,
                        sources={name: dict(sha256=checksum(path)) for name, path in
                                 [('mutations', mutations), ('classification', classification), ('vcf', vcf)]})
        db.execute('INSERT INTO metadata VALUES (?)', (json.dumps(metadata),))
        db.commit()
    except BaseException:
        db.close()
        temporary.unlink(missing_ok=True)
        raise
    db.close()
    temporary.rename(output)
    return metadata


class CosmicEvidence:
    def __init__(self, path):
        self.db = sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)
        self.db.row_factory = sqlite3.Row
        try:
            self.metadata = json.loads(self.db.execute('SELECT json FROM metadata').fetchone()[0])
            if self.metadata.get('schema_version') != 1 or self.metadata.get('assembly') != 'GRCh38':
                raise ValueError('Unsupported COSMIC evidence database')
            self.db.execute('SELECT ' + ','.join(EVIDENCE_FIELDS) + ' FROM evidence LIMIT 1').fetchone()
            self.db.execute('SELECT id,json FROM phenotypes LIMIT 1').fetchone()
        except Exception:
            self.db.close()
            raise

    def lookup(self, identifier):
        records = []
        phenotypes = {}
        samples, blood_samples = set(), set()
        for row in self.db.execute('SELECT e.*, p.json AS phenotype FROM evidence e LEFT JOIN phenotypes p ON e.COSMIC_PHENOTYPE_ID=p.id WHERE GENOMIC_MUTATION_ID=?', (identifier,)):
            evidence = {k: row[k] for k in EVIDENCE_FIELDS}
            phenotype = json.loads(row['phenotype']) if row['phenotype'] else None
            sample = evidence['COSMIC_SAMPLE_ID']
            valid_sample = sample not in ('', '.', 'NS')
            if valid_sample:
                samples.add(sample)
            if phenotype and phenotype['PRIMARY_SITE'] == BLOOD_SITE and valid_sample:
                blood_samples.add(sample)
            pid = evidence['COSMIC_PHENOTYPE_ID']
            if pid not in phenotypes:
                phenotypes[pid] = dict(classification=phenotype, samples=set(), publications=set(), somatic_statuses=set())
            group = phenotypes[pid]
            if valid_sample:
                group['samples'].add(sample)
            if evidence['PUBMED_PMID'] not in ('', '.', 'NS'):
                group['publications'].add(evidence['PUBMED_PMID'])
            group['somatic_statuses'].add(evidence['MUTATION_SOMATIC_STATUS'])
            records.append(evidence)
        return dict(identifier=identifier, status='matched' if records else 'missing_targeted_evidence',
                    sample_count=len(samples), blood_sample_count=len(blood_samples),
                    blood_match=any(g['classification'] and g['classification']['PRIMARY_SITE'] == BLOOD_SITE for g in phenotypes.values()),
                    phenotypes=[dict(phenotype_id=pid, classification=g['classification'],
                                     sample_count=len(g['samples']), publications=sorted(g['publications']),
                                     somatic_statuses=sorted(g['somatic_statuses']),
                                     status='matched' if g['classification'] else 'missing_classification')
                                for pid, g in sorted(phenotypes.items())], evidence=records)

    def close(self):
        self.db.close()
