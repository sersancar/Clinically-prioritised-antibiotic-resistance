#!/usr/bin/env python3
"""Pool qualifying counts over ALL three filters' post-QC read-end totals."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from ena_manifest import read_tsv,validate_filters


def pool(filters, measurements, totals):
    groups=validate_filters(filters,require_runs=False)
    metadata={r['filter_id']:r for r in filters}
    if set(totals)!=set(metadata) or any(v<=0 for v in totals.values()):
        raise ValueError('Require positive runtime post-QC totals for every filter')
    seen=set();data=defaultdict(list);references={}
    for row in measurements:
        fid=row['filter_id'];ref=row['reference'];key=(fid,ref)
        if fid not in metadata or key in seen:raise ValueError('Unknown/duplicate filter-reference')
        seen.add(key)
        length=int(row['length'])
        if ref in references and references[ref]!=length:raise ValueError('Inconsistent reference lengths')
        references[ref]=length
        data[(metadata[fid]['biosample'],ref)].append(row)
    # Metrics include every FASTA reference per filter, including zero alignments.
    if len(seen)!=len(filters)*len(references):raise ValueError('Incomplete filter-reference measurements')
    result=[]
    for bio,group in sorted(groups.items()):
        denominator=sum(totals[r['filter_id']] for r in group)
        for ref,length in references.items():
            rows=data[(bio,ref)]
            for screen in ['primary','relaxed']:
                qualifying=[r for r in rows if str(r[screen]).lower()=='true']
                count=sum(int(r['mapped_read_ends']) for r in qualifying)
                result.append(dict(biosample=bio,site=group[0]['site'],date=group[0]['date'],
                                   reference=ref,length=length,screen=screen,
                                   positive_filters=len(qualifying),counts=count,total_read_ends=denominator,
                                   RPM=count*1e6/denominator,RPKM=count*1e9/denominator/length))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--filters',required=True);p.add_argument('--metrics',nargs='+',required=True)
    p.add_argument('--fastp',nargs='+',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();totals={}
    for path in a.fastp:
        fid=Path(path).name.removesuffix('.fastp.json')
        if fid in totals:raise ValueError('Duplicate QC JSON')
        totals[fid]=int(json.loads(Path(path).read_text())['summary']['after_filtering']['total_reads'])
    rows=pool(read_tsv(a.filters),[r for p in a.metrics for r in read_tsv(p)],totals)
    with open(a.output,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
