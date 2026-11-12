#!/usr/bin/env python3
"""Streaming BAM counts and coverage with explicit reference-length denominators."""
import argparse
import subprocess
from collections import Counter
import pysam
from common import fasta_lengths,write_tsv


def metrics(bam,fasta,filter_id,catalogue):
    lengths=fasta_lengths(fasta);counts=Counter();covered=Counter();depthsum=Counter()
    with pysam.AlignmentFile(bam,'rb') as f:
        if dict(zip(f.references,f.lengths))!=lengths:raise ValueError('BAM and reference FASTA differ')
        for read in f.fetch(until_eof=True):
            if not read.flag & 2308 and read.mapping_quality>=30:counts[read.reference_name]+=1
    args=['samtools','mpileup','-aa','-B','-d','0','-q','30','-Q','0',bam] if catalogue=='crass' else \
         ['samtools','depth','-aa','-d','0','-q','0','-Q','30',bam]
    with subprocess.Popen(args,stdout=subprocess.PIPE,text=True) as process:
        for line in process.stdout:
            fields=line.rstrip().split('\t');ref=fields[0]
            if ref not in lengths:raise ValueError('Unknown coverage reference')
            value=int(fields[3] if catalogue=='crass' else fields[2])
            covered[ref]+=value>0;depthsum[ref]+=value
        if process.wait()!=0:raise RuntimeError('SAMtools coverage failed')
    return [dict(filter_id=filter_id,catalogue=catalogue,reference=ref,length=L,
                 covered_bases=covered[ref],breadth=covered[ref]/L,mean_depth=depthsum[ref]/L,
                 mapped_read_ends=counts[ref],primary=covered[ref]/L>=.8 and (catalogue=='broad' or depthsum[ref]/L>=5),
                 relaxed=covered[ref]/L>=.8) for ref,L in lengths.items()]


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['bam','fasta','filter-id','catalogue','output']:p.add_argument('--'+key,required=True)
    a=p.parse_args();write_tsv(a.output,metrics(a.bam,a.fasta,a.filter_id,a.catalogue))
