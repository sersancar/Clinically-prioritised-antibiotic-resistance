#!/usr/bin/env python3
"""Explicit, revision-checked launcher for the supplied AMR++ command.

Run this external Nextflow stage separately; do not hide 60 child workflows
inside the parent workflow. Its own work directory supports -resume.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


def build_command(executable,source,config,reads,output,profile,threads,work):
    if threads<1:raise ValueError('Threads must be positive')
    if not profile or not re.fullmatch('[A-Za-z0-9_,.-]+',profile):raise ValueError('Invalid profile')
    if not Path(reads).is_absolute():raise ValueError('Supply an absolute read glob')
    return [executable,'run',str(Path(source)/'main_AMR++.nf'),
            '-c',str(config),'-profile',profile,'-resume','-work-dir',str(work),
            '--snp','Y','--deduped','Y','--threads',str(threads),
            '--pipeline','resistome','--reads',reads,'--output',str(output)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','revision','config','reads-glob','output','profile','run-dir']:
        p.add_argument('--'+key,required=True)
    p.add_argument('--nextflow',default='nextflow');p.add_argument('--threads',type=int,default=36)
    p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    source=Path(a.source).resolve();config=Path(a.config).resolve();output=Path(a.output).resolve()
    run=Path(a.run_dir).resolve();run.mkdir(parents=True,exist_ok=True)
    if not re.fullmatch('[a-f0-9]{40}',a.revision):raise ValueError('Pin the full AMR++ Git commit')
    if not (source/'main_AMR++.nf').is_file() or not config.is_file():raise ValueError('AMR++ source/config missing')
    observed=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    if observed!=a.revision:raise ValueError('AMR++ revision mismatch')
    for args in [['diff','--quiet'],['diff','--cached','--quiet']]:
        subprocess.run(['git','-C',str(source),*args],check=True)
    command=build_command(a.nextflow,source,config,a.reads_glob,output,a.profile,a.threads,run/'work')
    configs={str(x.relative_to(source)):hashlib.sha256(x.read_bytes()).hexdigest()
             for x in sorted(source.rglob('*.config')) if '.git' not in x.parts and 'work' not in x.parts}
    record={'command':command,'source_revision':observed,'additional_config_sha256':hashlib.sha256(config.read_bytes()).hexdigest(),
            'source_configs_sha256':configs,'dry_run':a.dry_run,
            'resources':'Record reference FASTA/annotation/SNP checksums separately; configs alone are insufficient.',
            'source_command':'AMRplusplus_launch_ncahill.sh; --snp Y --deduped Y --pipeline resistome',
            'java':'Use the validated execution environment; do not silently replace its Java/Nextflow.'}
    (run/'amrplusplus_launch.json').write_text(json.dumps(record,indent=2)+'\n')
    if a.dry_run:
        print('Validated launch recorded; AMR++ was not executed.')
        return
    env=os.environ.copy();env['NXF_CLEANUP']='false'
    subprocess.run(command,cwd=run,env=env,check=True)


if __name__=='__main__':main()
