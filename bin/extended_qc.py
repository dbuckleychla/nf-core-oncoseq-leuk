#!/usr/bin/env python3
"""Additional sequencing QC without modifying native analyses."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qc_extended.contracts import FILES, checksum, identifier, state, write_json, write_tables


def load(path):
    return json.loads(Path(path).read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['reference', 'alignment', 'methylation', 'assemble', 'preprocessing', 'report', 'attach-report', 'run-summary'])
    p.add_argument('--fragment')
    for key in ('bam', 'fasta', 'reference', 'cpgs', 'enrichment', 'targets', 'alignment', 'methylation', 'bedmethyl', 'input', 'preprocessing', 'panel_manifest'):
        p.add_argument('--' + key.replace('_', '-'))
    p.add_argument('--sample')
    p.add_argument('--assembly', default='hg38')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--disable-methylation', action='store_true')
    p.add_argument('--code-digest')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    if a.threads < 1:
        p.error('threads must be positive')
    if a.command == 'reference':
        from qc_extended.regions import reference
        reference(a.fasta, a.output, a.cpgs)
    elif a.command == 'alignment':
        from qc_extended.alignment import collect
        data = collect(a.bam, load(a.reference), a.enrichment, a.targets, a.threads)
        mods = data['modifications']
        reason = ('user_disabled' if a.disable_methylation else
                  'invalid_modification_encoding' if mods.get('mapped_invalid') else
                  'unsupported_cytosine_modification' if mods.get('mapped_unsupported') else
                  'no_usable_modifications' if not mods.get('usable') else None)
        data['methylation_status'] = state('disabled' if a.disable_methylation else 'unavailable' if reason else 'partial' if mods.get('usable', 0) < data['counts'].get('mapped_primary', 0) else 'complete', reason or ('partial_modification_evidence' if mods.get('usable', 0) < data['counts'].get('mapped_primary', 0) else None))
        # Debian build-option text can contain non-UTF8 bytes; only the version
        # header is part of QC provenance, not compiler paths/options.
        version = subprocess.check_output(['samtools', '--version']).splitlines()[0].decode('ascii')
        data['versions'] = dict(samtools=version)
        write_json(a.output, data)
    elif a.command == 'methylation':
        from qc_extended.methylation import collect, panels
        aln, ref = load(a.alignment), load(a.reference)
        panel_data = panels(a.panel_manifest, a.fasta, ref['lengths'], a.assembly)
        status = aln['methylation_status']
        if status['status'] in ('disabled', 'unavailable'):
            features = [dict(**{k: v for k, v in panel.items() if k != 'features'}, size=len(panel['features']), depth=None)
                        for panel in panel_data]
            write_json(a.output, dict(status=status, scopes={}, features=features))
        else:
            supported = [code for code in ('m', 'h') if aln['modifications'].get('mapped_' + code)]
            result = collect(a.bedmethyl, ref, a.cpgs, a.enrichment, a.targets, supported, panel_data)
            result['status'] = status
            write_json(a.output, result)
    elif a.command == 'preprocessing':
        from qc_extended.preprocessing import collect
        write_json(a.output, collect(a.input))
    elif a.command == 'assemble':
        from qc_extended.assemble import assemble
        if a.preprocessing:
            stages = load(a.preprocessing)
        else:
            stages = []
        data = assemble(a.sample, load(a.alignment), load(a.methylation), load(a.reference),
                        stages,
                        dict(assembly=a.assembly, enrichment_sha256=checksum(a.enrichment), targets_sha256=checksum(a.targets),
                             code_sha256=a.code_digest,
                             images=dict(metrics='ghcr.io/dbuckleychla/frankontstein@sha256:15c315c121d37b5c115d3d965836e60d56756498ca2deebe066a373d41e96785',
                                         methylation='ghcr.io/chusj-pigu/classy-slim@sha256:7ec5fc1e4607be13d30dc65a1b07f2844f8289edbec0ca67a70e9b88ef8beb4c'),
                             versions=load(a.alignment).get('versions', {})))
        out = Path(a.output)
        out.mkdir(parents=True, exist_ok=True)
        write_json(out / 'metrics.json', data)
        write_tables(out, data)
    elif a.command == 'report':
        from qc_extended.report import content
        Path(a.output).write_text(content(load(a.input)))
    elif a.command == 'attach-report':
        import re
        text = Path(a.input).read_text()
        fragment = Path(a.fragment).read_text()
        marker = '<!-- ONCOSEQ_EXTENDED_QC -->'
        if marker in text:
            raise ValueError('Extended QC already attached to report')
        content = marker + '<section id="extended-qc">' + fragment + '</section>'
        if '</main>' not in text:
            raise ValueError('Native report has no main content element')
        text = text.replace('</main>', content + '\n</main>', 1)
        # Add navigation using the same Quarto TOC convention, if a TOC exists.
        navigation = re.search(r'(<nav[^>]*id="TOC"[\s\S]*?)(</ul>\s*</nav>)', text)
        if navigation:
            entry = '<li><a href="#extended-qc" class="nav-link">Extended sequencing QC</a></li>'
            text = text[:navigation.end(1)] + entry + text[navigation.end(1):]
        Path(a.output).write_text(text)
    else:
        import csv
        paths = sorted(Path(a.input).glob('*/summary.tsv'))
        with open(a.output, 'w') as out:
            writer = None
            for path in paths:
                with path.open() as source:
                    reader = csv.DictReader(source, delimiter='\t')
                    if writer is None:
                        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter='\t', lineterminator='\n')
                        writer.writeheader()
                    writer.writerows(reader)


if __name__ == '__main__':
    main()
