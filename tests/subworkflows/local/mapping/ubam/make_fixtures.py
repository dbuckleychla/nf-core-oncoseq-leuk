"""Small synthetic ONT-like uBAMs; no patient reads."""
import array
import random
from pathlib import Path
import pysam


def make(root):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    random.seed(1729)
    sequence=''.join(random.choice('ACGT') for _ in range(8000))
    (root/'reference.fa').write_text('>chr1\n'+sequence+'\n')
    (root/'reads.fastq').write_text('@SYNTHETIC_FASTQ\n'+sequence[200:1200]+'\n+\n'+'I'*1000+'\n')
    header={'HD':{'VN':'1.6','SO':'unknown'},'RG':[{'ID':'synthetic','SM':'SYNTHETIC'}]}
    for name in ('valid','missing_mm','missing_st','missing_after_sample'):
        with pysam.AlignmentFile(root/(name+'.bam'),'wb',header=header) as output:
            for i in range(1005):
                read=pysam.AlignedSegment()
                read.query_name=f'SYNTHETIC_{i:04d}';read.flag=4
                read.query_sequence=sequence[(i%25)*70:(i%25)*70+1000]
                read.query_qualities=pysam.qualitystring_to_array('I'*1000)
                if not (name=='missing_mm' and i==500) and not (name=='missing_after_sample' and i>=1000):
                    read.set_tag('MM','C+m,0;')
                read.set_tag('ML',array.array('B',[220]))
                read.set_tag('MN',1000)
                read.set_tag('qs',5.0 if i==0 else 20.0)
                read.set_tag('RG','synthetic')
                if name!='missing_st':read.set_tag('st','2020-01-01T00:00:00Z')
                output.write(read)


if __name__=='__main__':
    make(Path(__file__).resolve().parents[5]/'test_data/ubam')
