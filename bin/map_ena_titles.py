#!/usr/bin/env python3
"""Match explicit ENA sample titles to the manuscript's canonical filter IDs.

Freshwater/FW denotes site D; Seawater/SW prefixes retain their lettered site.
No numeric-alias or row-order inference is permitted.
"""
import argparse
import csv
import re
from ena_manifest import read_tsv,validate_filters


def canonical(title):
    match=re.fullmatch(r'([ABC])([123])(?:Seawater|SW)(\d{6})',title)
    if match:site,rep,date=match.groups()
    else:
        match=re.fullmatch(r'([123])(?:Freshwater|FW)(\d{6})',title)
        if not match:raise ValueError(f'Unrecognised ENA title: {title}')
        rep,date=match.groups();site='D'
    return f'{site}_{rep}_{date}'


def map_titles(filters,ena):
    mapped={}
    for row in ena:
        fid=canonical(row['sample_title'])
        if fid in mapped:raise ValueError(f'Multiple ENA runs for {fid}; resolve manually')
        mapped[fid]=row['run_accession']
    if set(mapped)!={r['filter_id'] for r in filters}:
        raise ValueError('ENA title set does not match the 60 manuscript filter IDs')
    output=[dict(r,run_accession=mapped[r['filter_id']]) for r in filters]
    validate_filters(output)
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ['filters','ena-report','output']:p.add_argument('--'+x,required=True)
    a=p.parse_args();rows=map_titles(read_tsv(a.filters),read_tsv(a.ena_report))
    with open(a.output,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
