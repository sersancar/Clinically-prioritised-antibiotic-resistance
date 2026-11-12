#!/usr/bin/env python3
"""Small deterministic paired reads, AMR references and taxonomy for tool tests."""
import gzip
import random
import sys
from pathlib import Path


def make(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);tax=out/'taxonomy';tax.mkdir(exist_ok=True)
    ranks=[(1,1,'no rank','root'),(2,1,'superkingdom','Bacteria'),(10,2,'class','Synthetic class'),
           (20,10,'genus','Synthetic genus'),(30,20,'species','Synthetic species A'),(31,20,'species','Synthetic species B')]
    (tax/'nodes.dmp').write_text(''.join(f'{i}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t11\t|\t1\t|\t1\t|\t0\t|\t0\t|\t0\t|\t\t|\n' for i,parent,rank,name in ranks))
    (tax/'names.dmp').write_text(''.join(f'{i}\t|\t{name}\t|\t\t|\tscientific name\t|\n' for i,parent,rank,name in ranks))
    random.seed(8723);sequences=[''.join(random.choice('ACGT') for _ in range(1200)) for _ in range(2)]
    complement=str.maketrans('ACGT','TGCA')
    (out/'genomes.fa').write_text(''.join(f'>synthetic{i}|kraken:taxid|{30+i}\n{seq}\n' for i,seq in enumerate(sequences)))
    (out/'megares.fa').write_text(''.join(f'>MEG_{i+1}|Drugs|MLS|Mechanism|gene{i+1}\n{seq}\n' for i,seq in enumerate(sequences)))
    (out/'carpdm.fa').write_text(''.join(f'>gb|synthetic{i}|+|1-1200|ARO:{i+1}|'+('sul2' if i==0 else 'tet(O)')+'\n'+seq+'\n' for i,seq in enumerate(sequences)))
    with gzip.open(out/'reads_1.fq.gz','wt') as r1,gzip.open(out/'reads_2.fq.gz','wt') as r2:
        for i,seq in enumerate(sequences):
            for pos in range(0,901,15):
                name=f'pair{i}_{pos}';a=seq[pos:pos+150];b=seq[pos+150:pos+300].translate(complement)[::-1]
                r1.write(f'@{name}/1\n{a}\n+\n'+('I'*150)+'\n');r2.write(f'@{name}/2\n{b}\n+\n'+('I'*150)+'\n')


if __name__=='__main__':make(sys.argv[1])
