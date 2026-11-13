#!/usr/bin/env python3
"""Prepare versioned reference resources before any sequencing download."""
import argparse
import gzip
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path
from common import fasta_lengths,sha256

AMR_COMMIT='fded6074d493fcdaf1ffd0c192b6b794249a464d'
SNP_COMMIT='8401ce3fff2860b257a9bcef23ac155e1e3e5222'


def retrieve(url,destination):
    with urllib.request.urlopen(url,timeout=120) as response,open(destination,'wb') as f:
        shutil.copyfileobj(response,f)


def copy_fasta(source,destination):
    """Write plain FASTA, detecting gzip even after Nextflow renames the input."""
    with open(source,'rb') as probe:
        compressed=probe.read(2)==b'\x1f\x8b'
    opener=gzip.open if compressed else open
    with opener(source,'rb') as src,open(destination,'wb') as dst:
        shutil.copyfileobj(src,dst)


def prepare(carpdm,megares,crass,snp_bundle,taxonomy,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);sources={}
    for name,source in [('carpdm',carpdm),('megares',megares),('crass',crass)]:
        dest=out/(name+'.fa')
        if source and Path(source).name!='NO_INPUT':
            copy_fasta(source,dest);sources[name]=str(Path(source).resolve())
        elif name=='megares':
            url=f'https://raw.githubusercontent.com/Microbial-Ecology-Group/AMRplusplus/{AMR_COMMIT}/data/amr/megares_database_v3.00.fasta'
            retrieve(url,dest);sources[name]=url
        elif name=='crass':
            url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_024711.1&rettype=fasta&retmode=text'
            retrieve(url,dest);sources[name]=url
        else:raise ValueError('Supply the CARPDM nucleotide reference FASTA')
        lengths=fasta_lengths(dest)
        if name=='crass' and (len(lengths)!=1 or not any('NC_024711.1' in x for x in lengths)):
            raise ValueError('Marker FASTA must contain only NC_024711.1')
    if snp_bundle and Path(snp_bundle).name!='NO_INPUT':
        shutil.copytree(snp_bundle,out/'snp',symlinks=False,ignore=shutil.ignore_patterns('.git','__pycache__','*.pyc'))
        sources['snp']=str(Path(snp_bundle).resolve())
    else:
        subprocess.run(['git','init',str(out/'snp')],check=True,capture_output=True)
        subprocess.run(['git','-C',str(out/'snp'),'fetch','--depth','1',
                        'https://github.com/Microbial-Ecology-Group/AmrPlusPlus_SNP.git',SNP_COMMIT],check=True)
        subprocess.run(['git','-C',str(out/'snp'),'checkout','--detach','FETCH_HEAD'],check=True)
        shutil.rmtree(out/'snp'/'.git');sources['snp']=SNP_COMMIT
    for required in ['SNP_Verification.py','config.ini','data/SNPinfo.fasta']:
        if not (out/'snp'/required).is_file():raise ValueError(f'SNP bundle missing {required}')
    tax=Path(taxonomy)
    for name in ['nodes.dmp','names.dmp']:
        if not (tax/name).is_file():raise ValueError(f'Frozen taxonomy missing {name}')
    files={p.relative_to(out).as_posix():sha256(p) for p in out.rglob('*') if p.is_file()}
    files.update({f'taxonomy/{name}':sha256(tax/name) for name in ['nodes.dmp','names.dmp']})
    (out/'provenance.json').write_text(json.dumps({'sources':sources,'sha256':files,
        'megares_default_version':'3.00','snp_default_commit':SNP_COMMIT},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['carpdm','megares','crass','snp-bundle','taxonomy','out']:p.add_argument('--'+key,required=True)
    a=p.parse_args();prepare(a.carpdm,a.megares,a.crass,a.snp_bundle,a.taxonomy,a.out)
