#!/usr/bin/env python3
"""Resolve explicitly mapped ENA runs without inferring sample aliases."""
import argparse
import csv
import io
import json
import re
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

FIELDS = 'run_accession,sample_accession,sample_alias,sample_title,library_layout,fastq_ftp,fastq_md5,fastq_bytes'


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def validate_filters(rows, require_runs=True):
    if len(rows) != 60:
        raise ValueError('Paper mode requires exactly 60 filters')
    ids = [r['filter_id'] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate filter_id')
    groups = {}
    for r in rows:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', r['filter_id']):
            raise ValueError('Unsafe filter identifier')
        if r['site'] not in 'ABCD' or len(r['site']) != 1:
            raise ValueError('Site must be A, B, C or D')
        if r['biosample'] != f"{r['site']}_{r['date']}":
            raise ValueError('Biological sample does not match site/date')
        group = groups.setdefault(r['biosample'], [])
        group.append(r)
        runs = r['run_accession'].split(';')
        if require_runs and not all(re.fullmatch(r'(ERR|SRR|DRR)\d+', x) for x in runs):
            raise ValueError(f"Verified ENA run mapping missing/invalid: {r['filter_id']}")
    if len(groups) != 20:
        raise ValueError('Expected 20 biological samples')
    for key, group in groups.items():
        if sorted(r['replicate'] for r in group) != ['1', '2', '3']:
            raise ValueError(f'{key}: expected filters 1, 2, 3')
        for field in ['salinity', 'pH', 'water_temperature']:
            if len({r[field] for r in group}) != 1:
                raise ValueError(f'{key}: inconsistent biological metadata {field}')
    expected = {'2021-03-29','2021-03-31','2021-04-12','2021-04-14','2021-12-02','2021-12-06'}
    for site in 'ABCD':
        dates = {r['date'] for r in rows if r['site'] == site}
        if dates != (expected & {'2021-12-02','2021-12-06'} if site == 'C' else expected):
            raise ValueError(f'Unexpected sampling dates at {site}')
    all_runs = [x for r in rows for x in r['run_accession'].split(';') if x]
    if len(set(all_runs)) != len(all_runs):
        raise ValueError('Run assigned to multiple filters')
    return groups


def resolve(filters, ena):
    validate_filters(filters)
    indexed = {r['run_accession']: r for r in ena}
    if len(indexed) != len(ena):
        raise ValueError('Duplicate ENA run records')
    result = []
    for row in filters:
        for run in row['run_accession'].split(';'):
            if run not in indexed:
                raise ValueError(f'{run} not present in requested study')
            record = indexed[run]
            urls = record['fastq_ftp'].split(';')
            md5s = record['fastq_md5'].split(';')
            if record['library_layout'] != 'PAIRED' or len(urls) != 2 or len(md5s) != 2:
                raise ValueError(f'{run}: expected exactly two paired FASTQ files')
            for url, md5 in zip(urls, md5s):
                parsed = urllib.parse.urlsplit('https://' + url.removeprefix('ftp://').removeprefix('https://'))
                if parsed.hostname != 'ftp.sra.ebi.ac.uk' or not re.fullmatch('[0-9a-fA-F]{32}', md5):
                    raise ValueError(f'{run}: invalid ENA URL/checksum')
            if not urls[0].endswith('_1.fastq.gz') or not urls[1].endswith('_2.fastq.gz'):
                raise ValueError(f'{run}: mate order needs manual resolution')
            result.append(dict(filter_id=row['filter_id'],run_accession=run,
                               url1='https://'+urls[0].removeprefix('ftp://').removeprefix('https://'),
                               url2='https://'+urls[1].removeprefix('ftp://').removeprefix('https://'),
                               md5_1=md5s[0],md5_2=md5s[1]))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--filters', required=True)
    p.add_argument('--study', default='PRJEB75359')
    p.add_argument('--ena-report', help='Use a saved ENA TSV instead of network access')
    p.add_argument('--outdir', default='ena')
    a = p.parse_args()
    filters = read_tsv(a.filters)
    validate_filters(filters)  # fail before network activity when mapping is incomplete
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    if a.ena_report:
        raw = Path(a.ena_report).read_text()
    else:
        query = urllib.parse.urlencode(dict(accession=a.study,result='read_run',fields=FIELDS,format='tsv'))
        with urllib.request.urlopen('https://www.ebi.ac.uk/ena/portal/api/filereport?'+query,timeout=60) as r:
            raw = r.read().decode()
    records = resolve(filters,list(csv.DictReader(io.StringIO(raw),delimiter='\t')))
    (out/'ena_report.tsv').write_text(raw)
    with (out/'downloads.tsv').open('w') as f:
        writer = csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t')
        writer.writeheader();writer.writerows(records)
    (out/'selection.json').write_text(json.dumps({'study':a.study,'filters':len(filters),'runs':len(records)},indent=2)+'\n')


if __name__ == '__main__':
    main()
