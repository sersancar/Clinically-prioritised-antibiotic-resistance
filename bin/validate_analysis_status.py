#!/usr/bin/env python3
"""Require the critical original R modules; completion alone is insufficient."""
import csv
import sys

if __name__=='__main__':
    with open(sys.argv[1]) as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
    status={r['Module']:r['Status'] for r in rows}
    required={'replication':'DONE','length_normalization':'DONE','CARPDM':'DONE',
              'taxonomy':'DESCRIPTIVE','environment':'DONE','pipeline':'COMPLETE'}
    bad=[f'{k}: {status.get(k,"missing")}' for k,v in required.items() if status.get(k)!=v]
    if bad:raise SystemExit('Incomplete original R analysis: '+', '.join(bad))
