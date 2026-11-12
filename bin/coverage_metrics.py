#!/usr/bin/env python3
"""Summarise SAMtools depth output using full FASTA reference lengths.

CARPDM input: samtools depth -aa -q 0 -Q 30 (no -s; no overlap removal).
Marker input: first four fields of samtools mpileup -aa -B -d 0 -q 30 -Q 0.
The SAM input for counts must have the same MAPQ/flag filtering as the BAM.
"""
import argparse
import csv
from collections import Counter


def fasta_lengths(path):
    lengths={};name=None
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in lengths: raise ValueError('Duplicate FASTA reference')
                lengths[name]=0
            elif line.strip():
                if name is None: raise ValueError('Sequence before FASTA header')
                lengths[name]+=len(line.strip())
    if not lengths or min(lengths.values())<=0: raise ValueError('Empty reference sequence')
    return lengths


def summarise(fasta, depth, sam, pileup=False):
    lengths=fasta_lengths(fasta);covered=Counter();total=Counter();counts=Counter()
    # Depth is streamed; no genome-sized arrays and no BAM-sized identifier set.
    previous_ref=None; previous_pos=0; finished=set()
    with open(depth) as f:
        for line in f:
            fields=line.rstrip().split('\t');ref=fields[0];pos=int(fields[1])
            value=int(fields[3 if pileup else 2])
            if ref not in lengths or not 1<=pos<=lengths[ref] or value<0:
                raise ValueError('Invalid depth record')
            if ref!=previous_ref:
                if ref in finished: raise ValueError('Depth reference appears in multiple blocks')
                if previous_ref is not None: finished.add(previous_ref)
                previous_ref=ref;previous_pos=0
            if pos<=previous_pos: raise ValueError('Duplicate/unsorted depth position')
            previous_pos=pos;covered[ref]+=value>0;total[ref]+=value
    with open(sam) as f:
        for line in f:
            if line.startswith('@'):continue
            fields=line.split('\t');flag=int(fields[1]);ref=fields[2]
            if flag & 2308 or int(fields[4])<30: raise ValueError('SAM was not correctly filtered')
            if ref not in lengths:raise ValueError('Alignment reference not in FASTA')
            counts[ref]+=1
    return [dict(reference=ref,length=length,covered_bases=covered[ref],
                 breadth=covered[ref]/length,mean_depth=total[ref]/length,
                 mapped_read_ends=counts[ref],primary=covered[ref]/length>=.8 and total[ref]/length>=5,
                 relaxed=covered[ref]/length>=.8)
            for ref,length in lengths.items()]


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ['fasta','depth','sam','filter-id','output']:p.add_argument('--'+x,required=True)
    p.add_argument('--pileup',action='store_true');a=p.parse_args()
    rows=summarise(a.fasta,a.depth,a.sam,a.pileup)
    with open(a.output,'w') as f:
        w=csv.DictWriter(f,fieldnames=['filter_id',*rows[0]],delimiter='\t',lineterminator='\n');w.writeheader()
        w.writerows(dict(filter_id=a.filter_id,**r) for r in rows)
