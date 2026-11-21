"""Shared GFF gene/transcript overlap for native genomic breakpoints."""
import gzip
from collections import Counter
from urllib.parse import unquote
from check_reference import contigs


def attributes(text):
    return {k:unquote(v) for k,v in (part.split('=',1) for part in text.split(';') if '=' in part)}


def gene_index(gff,fai,db):
    lengths=contigs(fai);aliases={c:c for c in lengths}
    for c in lengths:
        aliases[c.removeprefix('chr')]=c
        if c in ('chrM','MT','M'):
            for name in ('chrM','MT','M'):aliases[name]=c
    db.executescript('CREATE TABLE genes (id TEXT PRIMARY KEY, name TEXT, chrom TEXT, start INTEGER, end INTEGER, strand TEXT); CREATE INDEX gene_span ON genes(chrom,start,end); CREATE TABLE transcripts (id TEXT PRIMARY KEY,parent TEXT,chrom TEXT,start INTEGER,end INTEGER,strand TEXT); CREATE TABLE exons (parent TEXT,start INTEGER,end INTEGER); CREATE INDEX transcript_parent ON transcripts(parent); CREATE INDEX exon_parent ON exons(parent);')
    opener=gzip.open if str(gff).endswith('.gz') else open
    def features():
        with opener(gff,'rt') as stream:
            for line in stream:
                if line.startswith('#'):continue
                f=line.rstrip().split('\t')
                if len(f)==9:yield f,attributes(f[8])
    for f,a in features():
        if f[2]=='region' and a.get('chromosome') in aliases:aliases[f[0]]=aliases[a['chromosome']]
    excluded=Counter()
    for f,a in features():
        c=aliases.get(f[0])
        if c is None:excluded[f[0]]+=1;continue
        start,end=int(f[3]),int(f[4])
        if not 1<=start<=end<=lengths[c]:raise ValueError('GFF feature outside reference: '+f[0])
        if f[2]=='gene' and a.get('ID'):
            db.execute('INSERT INTO genes VALUES (?,?,?,?,?,?)',(a['ID'],a.get('Name') or a.get('gene') or a['ID'],c,start,end,f[6]))
        elif f[2] in ('mRNA','transcript','lnc_RNA','ncRNA','rRNA','tRNA','primary_transcript') and a.get('ID'):
            for parent in a.get('Parent','').split(','):db.execute('INSERT OR IGNORE INTO transcripts VALUES (?,?,?,?,?,?)',(a['ID'],parent,c,start,end,f[6]))
        elif f[2]=='exon':
            for parent in a.get('Parent','').split(','):db.execute('INSERT INTO exons VALUES (?,?,?)',(parent,start,end))
    db.commit()
    if not db.execute('SELECT COUNT(*) FROM genes').fetchone()[0]:raise ValueError('No GFF genes match the alignment reference')
    return dict(excluded_contigs=dict(excluded),contig_mapping=aliases,contig_lengths=lengths)


def overlap(db,chrom,pos):
    if chrom is None or pos is None:return []
    result=[]
    for gid,name,start,end,strand in db.execute('SELECT id,name,start,end,strand FROM genes WHERE chrom=? AND start<=? AND end>=? ORDER BY name,id',(chrom,pos,pos)):
        transcripts=[]
        for tid,ts,te,tstrand in db.execute('SELECT id,start,end,strand FROM transcripts WHERE parent=? AND start<=? AND end>=? ORDER BY id',(gid,pos,pos)):
            exon=bool(db.execute('SELECT 1 FROM exons WHERE parent=? AND start<=? AND end>=? LIMIT 1',(tid,pos,pos)).fetchone())
            has_exons=bool(db.execute('SELECT 1 FROM exons WHERE parent=? LIMIT 1',(tid,)).fetchone())
            transcripts.append(dict(id=tid,strand=tstrand,context='exon' if exon else 'intron' if has_exons else 'transcript_overlap'))
        result.append(dict(id=gid,name=name,strand=strand,transcripts=transcripts,context='transcribed' if transcripts else 'gene_overlap'))
    return result
