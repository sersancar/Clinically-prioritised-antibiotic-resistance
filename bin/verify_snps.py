#!/usr/bin/env python3
"""Apply the pinned AMR++ SNP verifier to eligible, deduplicated read-end counts."""
import argparse
import configparser
import csv
import json
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO
from common import read_tsv,write_tsv


def verify(bam,metrics,bundle,fasta,filter_id,threads,out):
    rows=read_tsv(metrics);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    baseline={r['reference']:int(r['mapped_read_ends']) if r['primary']=='True' else 0 for r in rows}
    needed=any(n>0 and 'RequiresSNPConfirmation' in r for r,n in baseline.items())
    counts=baseline.copy()
    if needed:
        bundle=Path(bundle).resolve();fasta=Path(fasta).resolve()
        records=[];covered=set()
        with (bundle/'data/SNPinfo.fasta').open() as handle:
            for record in SeqIO.parse(handle,'fasta'):
                key='|'.join(record.id.split('|')[:5])+'|RequiresSNPConfirmation'
                if baseline.get(key,0)>0:records.append(record);covered.add(key)
        for ref in baseline:
            if 'RequiresSNPConfirmation' in ref and ref not in covered:counts[ref]=0
        needed=bool(records)
        snp_info=(out/'eligible_SNPinfo.fasta').absolute()
        SeqIO.write(records,snp_info,'fasta')
    if needed:
        local_bam=(out/(filter_id+'.bam')).absolute();local_bam.unlink(missing_ok=True)
        local_bam.symlink_to(Path(bam).resolve())
        matrix=(out/'input.csv').absolute()
        with matrix.open('w') as f:
            w=csv.writer(f);w.writerow(['gene_accession',filter_id]);w.writerows(baseline.items())
        config=configparser.ConfigParser();config.read(bundle/'config.ini')
        config['SOURCE_FILES']['megares_fasta']=str(fasta)
        config['SOURCE_FILES']['snp_info_fasta']=str(snp_info)
        config['SOURCE_FILES']['count_matrix']=str(matrix)
        config['SOURCE_FILES']['bam_input']=str(local_bam)
        config['FULL_FILE_NAMES']['sample']=filter_id
        config['SETTINGS']['reads']='false';config['SETTINGS']['detailed']='false'
        config['SETTINGS']['amrplusplus']='true';config['SETTINGS']['threads']=str(threads)
        config['FOLDERS']['main_output_folder']=str((out/'verified').absolute())+'/'
        config_path=(out/'verification.ini').absolute()
        with config_path.open('w') as f:config.write(f)
        command=[sys.executable,str(Path(__file__).with_name('snp_worker.py')),str(bundle),'-c',str(config_path),
                 '-a','true','-i',str(local_bam),'-o',str((out/'verified').absolute())+'/',
                 '-t',str(threads),'--count_matrix',str(matrix),
                 '--count_matrix_final','verified_counts.csv','--output_reads=false','--detailed_output=false']
        subprocess.run(command,check=True)
        files=list((out/'verified').glob('*verified_counts.csv'))
        if len(files)!=1:raise ValueError('SNP verifier did not produce one confirmed matrix')
        with files[0].open() as f:
            result=list(csv.DictReader(f))
        if {r['gene_accession'] for r in result}!=set(baseline):raise ValueError('SNP matrix reference set changed')
        for r in result:
            value=float(r[filter_id]);ref=r['gene_accession']
            if not value.is_integer() or not 0<=value<=baseline[ref]:raise ValueError('Invalid SNP-confirmed count')
            counts[ref]=int(value)
    output=[]
    for r in rows:
        ref=r['reference'];output.append(dict(filter_id=filter_id,reference=ref,length=int(r['length']),
            breadth=float(r['breadth']),mean_depth=float(r['mean_depth']),raw_count=baseline[ref],
            count=counts[ref],requires_snp='RequiresSNPConfirmation' in ref))
    write_tsv(out/'counts.tsv',output)
    (out/'status.json').write_text(json.dumps({'filter_id':filter_id,'snp_verification_enabled':True,
       'verifier_executed':needed,'positive_snp_targets_before':sum(n>0 and 'RequiresSNPConfirmation' in r for r,n in baseline.items()),
       'positive_snp_targets_after':sum(n>0 and 'RequiresSNPConfirmation' in r for r,n in counts.items())},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['bam','metrics','bundle','fasta','filter-id','out']:p.add_argument('--'+key,required=True)
    p.add_argument('--threads',type=int,default=4);a=p.parse_args()
    verify(a.bam,a.metrics,a.bundle,a.fasta,a.filter_id,a.threads,a.out)
