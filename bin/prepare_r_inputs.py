#!/usr/bin/env python3
"""Translate validated pipeline products to the original R input contract.

Never reads published read totals as runtime denominators. Requires every filter.
AMR counts are already MAPQ/breadth filtered upstream; this script cannot add
missing SNP screening or deduplication to those counts.
"""
import argparse
import csv
import json
import shutil
from pathlib import Path
from coverage_metrics import fasta_lengths
from ena_manifest import read_tsv,validate_filters


def legacy_id(row):
    prefix='F' if row['site']=='D' else row['site']
    water='FW' if row['site']=='D' else 'SW'
    date=row['date'][8:10]+row['date'][5:7]+row['date'][2:4]
    return prefix+row['replicate']+water+date


def write_tsv(path,rows):
    with open(path,'w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def prepare(filters,amr_dir,metrics_dir,tax_dir,fastp_dir,fasta,outdir):
    validate_filters(filters,require_runs=False)
    out=Path(outdir)
    if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be empty')
    out.mkdir(parents=True,exist_ok=True)
    (out/'CARD_genes').mkdir();(out/'CARPDM_taxonomy').mkdir()
    lengths=fasta_lengths(fasta);counts={};all_refs=set();den=[];meta=[];manifest=[]
    for row in filters:
        fid=row['filter_id'];sid=legacy_id(row)
        qc=json.loads((Path(fastp_dir)/(fid+'.fastp.json')).read_text())
        total=int(qc['summary']['after_filtering']['total_reads'])
        if total<=0:raise ValueError(f'{fid}: invalid runtime total')
        den.append(dict(sample=sid,total_reads=total))
        manifest.append(dict(SampleID=sid,Site=row['site'],Date=row['date'],BioSample=row['biosample']))
        meta.append(dict(Sample_name=sid,Salinity=row['salinity'],pH=row['pH'],
                         Sample_temperature=row['water_temperature'],Temperature=row['air_temperature'],
                         Wind=row['wind'],Max_UV_index=row['uv'],Rel_humidity=row['humidity'],
                         Rainfall=row['rainfall'],Tidal_level=row['tide']))
        count_file=Path(amr_dir)/(fid+'.amr')/'reference_counts.tsv'
        values={}
        for r in read_tsv(count_file):
            ref=r['gene_accession'];n=int(r['count'])
            if ref not in lengths or ref in values or n<0:raise ValueError(f'{fid}: invalid AMR reference/count')
            values[ref]=n
        counts[sid]=values;all_refs.update(k for k,v in values.items() if v>0)
        metric_file=Path(metrics_dir)/(fid+'.carpdm.metrics.tsv')
        metrics=read_tsv(metric_file)
        if not metrics:raise ValueError(f'{fid}: missing all-reference CARPDM measurements')
        targets=[];seen=set()
        for r in metrics:
            ref=r['reference'];L=int(r['length']);covered=int(r['covered_bases']);depth=float(r['mean_depth']);n=int(r['mapped_read_ends'])
            if r['filter_id']!=fid or ref in seen or L<=0 or not 0<=covered<=L or depth<0 or n<0:
                raise ValueError(f'{fid}: invalid CARPDM metric')
            seen.add(ref)
            targets.append(dict(gene=ref,length=L,cov_bases=covered,breadth=float(r['breadth']),mean_depth=depth,num_reads=n))
        write_tsv(out/'CARD_genes'/(sid+'.ARG_metrics.tsv'),targets)
        tax=Path(tax_dir)/(fid+'.CARD_ARG_taxonomy.tsv')
        # Confirm schema even for header-only zero-call filters.
        with tax.open() as f:
            header=next(csv.reader(f,delimiter='\t'))
        if not {'ARG_ID','read_count','genus'}.issubset(header):raise ValueError(f'{fid}: invalid taxonomy schema')
        shutil.copyfile(tax,out/'CARPDM_taxonomy'/(sid+'.CARD_ARG_taxonomy.tsv'))
    if not all_refs:raise ValueError('No positive AMR references')
    write_tsv(out/'metadata.tsv',meta);write_tsv(out/'sample_manifest.tsv',manifest)
    write_tsv(out/'total_reads.tsv',den)
    write_tsv(out/'megares_lengths.tsv',[dict(gene_accession=k,length=lengths[k]) for k in sorted(all_refs)])
    with (out/'dedup_AMR_analytic_matrix2.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');ids=[legacy_id(r) for r in filters]
        writer.writerow(['gene_accession',*ids])
        writer.writerows([ref,*[counts[sid].get(ref,0) for sid in ids]] for ref in sorted(all_refs))
    (out/'PROVENANCE.json').write_text(json.dumps({'filters':60,'biological_samples':20,
       'denominator':'fastp after_filtering.total_reads from supplied runtime JSONs',
       'carpdm':'All-reference metrics required; set carpdm_unfiltered only after provenance validation',
       'amr':'Already-filtered source counts; SNP/dedup provenance remains an upstream requirement'},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['filters','amr-dir','metrics-dir','tax-dir','fastp-dir','megares-fasta','outdir']:
        p.add_argument('--'+name,required=True)
    a=p.parse_args()
    prepare(read_tsv(a.filters),a.amr_dir,a.metrics_dir,a.tax_dir,a.fastp_dir,a.megares_fasta,a.outdir)
