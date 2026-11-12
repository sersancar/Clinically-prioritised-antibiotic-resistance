#!/usr/bin/env python3
"""Rebuild reference counts under one explicitly selected original policy.

Input depth must be generated from the BAM using the matching source command.
The original scripts round breadth percentage to two decimals BEFORE screening.
No SNP gate is introduced here: input BAM/upstream SNP provenance is required.
"""
import argparse
import csv
import json
import sys
from contextlib import nullcontext
from collections import Counter
from pathlib import Path

POLICIES={'legacy2308':{'count_exclude':2308,'depth_add_exclude':2048},
          'exclude_duplicates3332':{'count_exclude':3332,'depth_add_exclude':3332}}


def count_references(depth,sam,policy):
    if policy not in POLICIES:raise ValueError('Select an explicit source policy')
    mask=POLICIES[policy]['count_exclude'];lengths=Counter();covered=Counter();counts=Counter()
    previous=None;last=0;finished=set()
    with open(depth) as f:
        for line in f:
            ref,pos,value=line.rstrip().split('\t');pos=int(pos);value=int(value)
            if ref!=previous:
                if ref in finished:raise ValueError('Unsorted depth references')
                if previous is not None:finished.add(previous)
                previous=ref;last=0
            # -aa is mandatory. Reject partial records rather than changing the
            # reference-length denominator as the original AWK otherwise would.
            if pos!=last+1 or value<0:raise ValueError('Require contiguous -aa depth records')
            last=pos;lengths[ref]+=1;covered[ref]+=value>0
    with (nullcontext(sys.stdin) if str(sam)=='-' else open(sam)) as f:
        for line in f:
            if line.startswith('@'):continue
            fields=line.split('\t');flag=int(fields[1])
            if flag & mask or int(fields[4])<30:continue
            ref=fields[2]
            if ref not in lengths:raise ValueError('SAM reference missing from depth output')
            counts[ref]+=1
    metrics=[dict(gene_accession=ref,length=L,covered_bases=covered[ref],
                  breadth_percent_rounded=float(f'{100*covered[ref]/L:.2f}'),
                  count=counts[ref]) for ref,L in lengths.items()]
    accepted=[{'gene_accession':r['gene_accession'],'count':r['count']} for r in metrics
              if r['count']>0 and r['breadth_percent_rounded']>=80]
    return metrics,accepted


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['depth','sam','output-dir']:p.add_argument('--'+key,required=True)
    p.add_argument('--policy',choices=list(POLICIES),required=True);a=p.parse_args()
    metrics,accepted=count_references(a.depth,a.sam,a.policy)
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    for name,rows,fields in [('reference_counts.tsv',accepted,['gene_accession','count']),
                             ('reference_metrics.tsv',metrics,['gene_accession','length','covered_bases','breadth_percent_rounded','count'])]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    (out/'provenance.json').write_text(json.dumps({'policy':a.policy,**POLICIES[a.policy],
       'mapq':30,'breadth_gate':'Two-decimal percent >=80, reproducing the original AWK order',
       'snp_gate':'None in the supplied BAM-count rebuilding scripts; verify upstream effects'},indent=2)+'\n')
