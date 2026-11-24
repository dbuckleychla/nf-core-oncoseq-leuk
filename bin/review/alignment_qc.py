"""Primary-read alignment metrics; never divide by supplementary records."""
import json
import re
from pathlib import Path


def primary_counts(mapped, unmapped, source, scope='Mapping-stage BAM before downstream filtering'):
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in (mapped, unmapped)):
        raise ValueError('Primary read counts must be non-negative integers')
    total = mapped + unmapped
    return dict(schema_version=1, primary_reads=total, mapped_primary_reads=mapped,
                unmapped_primary_reads=unmapped, alignment_rate_percent=100*mapped/total if total else None,
                source=source, scope=scope)


def load_alignment(path):
    data=json.loads(Path(path).read_text())
    if data.get('schema_version') == 1 and 'mapped_primary_reads' in data:
        return primary_counts(data['mapped_primary_reads'],data['unmapped_primary_reads'],data.get('source','Primary counts'),data.get('scope','Mapping-stage BAM'))
    total=mapped=0
    for category in ('QC-passed reads','QC-failed reads'):
        values=data.get(category,{})
        if 'primary' not in values or 'primary mapped' not in values:
            raise ValueError('Flagstat output lacks primary-read counts')
        total+=values['primary']; mapped+=values['primary mapped']
    if mapped>total:
        raise ValueError('Mapped primary reads exceed the primary-read total')
    return primary_counts(mapped,total-mapped,'samtools flagstat; QC-passed and QC-failed primary reads')


def apply_alignment(data, qc):
    data['alignment_qc']=qc
    names={'Read alignment rate (%)','Mapped primary reads','Unmapped primary reads','Primary reads submitted to alignment'}
    metrics=[dict(metric='Read alignment rate (%)',value=round(qc['alignment_rate_percent'],2) if qc['alignment_rate_percent'] is not None else 'Not available'),
             dict(metric='Mapped primary reads',value=qc['mapped_primary_reads']),
             dict(metric='Unmapped primary reads',value=qc['unmapped_primary_reads']),
             dict(metric='Primary reads submitted to alignment',value=qc['primary_reads'])]
    data['qc_summary']=metrics+[r for r in data.get('qc_summary',[]) if r['metric'] not in names]
    title='Read alignment (' + (f"{qc['alignment_rate_percent']:.2f}% mapped" if qc['alignment_rate_percent'] is not None else 'rate unavailable') + ')'
    plot=dict(title=title,unit='Primary reads',bars=[dict(label='Mapped',value=qc['mapped_primary_reads']),dict(label='Unmapped',value=qc['unmapped_primary_reads'])])
    data['qc_plots']=[plot]+[p for p in data.get('qc_plots',[]) if not p['title'].startswith('Read alignment')]


def patch_html(input, output, alignment_qc):
    """Update only embedded QC data; retain all existing figures and report text."""
    text=Path(input).read_text()
    match=re.search(r'(<script[^>]*id="review-data"[^>]*>)(.*?)(</script>)',text,re.S)
    if not match:
        raise ValueError('No review-data payload found')
    data=json.loads(match[2]);apply_alignment(data,load_alignment(alignment_qc))
    payload=json.dumps(data,separators=(',',':'),allow_nan=False).replace('<','\\u003c').replace('&','\\u0026')
    Path(output).write_text(text[:match.start(2)]+payload+text[match.end(2):])
