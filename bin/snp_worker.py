#!/usr/bin/env python3
"""Run pinned AMR++ SNP scoring with a serial, socket-free gene scheduler.

Only orchestration changes: parse_config, variant construction, read verification
and count adjustment are the upstream implementation. Nextflow parallelises
libraries, avoiding an additional Manager server and process pool per library.
"""
import importlib.util
import sys
import threading
from pathlib import Path
import pysam


def main():
    bundle=Path(sys.argv[1]).resolve();arguments=sys.argv[2:]
    sys.path.insert(0,str(bundle))
    spec=importlib.util.spec_from_file_location('coastal_amr_snp',bundle/'SNP_Verification.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def serial_genes(config,arg_list,genes):
        path=config['TEMP_FILES']['TEMP_BAM_SORTED']
        pysam.sort('-o',path,config['SOURCE_FILES']['BAM_INPUT']);pysam.index(path,path+'.bai')
        lock=threading.Lock()
        return [module.iterate((name,record,config,arg_list,lock)) for name,record in genes.items()]
    module.process_genes=serial_genes
    sys.argv=[str(bundle/'SNP_Verification.py'),*arguments]
    module.main()


if __name__=='__main__':main()
