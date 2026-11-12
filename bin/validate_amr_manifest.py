#!/usr/bin/env python3
"""Validate a complete canonical filter-to-deduplicated-BAM crosswalk."""
import argparse
import csv
from pathlib import Path
from ena_manifest import read_tsv,validate_filters


def validate(filters,rows,base):
    validate_filters(filters,require_runs=False)
    expected={r['filter_id'] for r in filters}
    if len(rows)!=60 or {r['filter_id'] for r in rows}!=expected:
        raise ValueError('Require one BAM for each of the 60 manuscript filters')
    paths=set();out=[]
    for r in rows:
        p=Path(r['bam']);p=(base/p).resolve() if not p.is_absolute() else p.resolve()
        if not p.is_file():raise ValueError(f'Missing BAM for {r["filter_id"]}: {p}')
        if p in paths:raise ValueError('Same BAM assigned to more than one filter')
        paths.add(p);out.append(dict(filter_id=r['filter_id'],bam=str(p)))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ['filters','manifest','output']:p.add_argument('--'+x,required=True)
    a=p.parse_args();manifest=Path(a.manifest).resolve()
    # Resolve relative BAM paths against the original manifest's location. When
    # invoked by Nextflow, the manifest is a symlink to that original file.
    rows=validate(read_tsv(a.filters),read_tsv(manifest),manifest.parent)
    with open(a.output,'w') as f:
        w=csv.DictWriter(f,fieldnames=['filter_id','bam'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
