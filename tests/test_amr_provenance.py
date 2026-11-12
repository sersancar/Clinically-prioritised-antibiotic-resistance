import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from count_amr_bam import count_references
from run_amrplusplus import build_command
from ena_manifest import read_tsv
from validate_amr_manifest import validate


class AMRProvenance(unittest.TestCase):
    def test_launcher_preserves_supplied_flags(self):
        command=build_command('nextflow','/amr','/config','/reads/*_clean_{1,2}.fq.gz','/out','local_slurm',36,'/work')
        for option,value in [('--snp','Y'),('--deduped','Y'),('--pipeline','resistome'),('--threads','36')]:
            self.assertEqual(command[command.index(option)+1],value)
        self.assertIn('-resume',command)
        with self.assertRaises(ValueError):build_command('nf','/amr','/conf','relative/*','/out','local',1,'/work')

    def test_duplicate_flags_change_counts_as_in_originals(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'depth').write_text(''.join(f'g\t{i}\t1\n' for i in range(1,11)))
            # One primary, one DUP, plus secondary/unmapped and low-MAPQ records.
            records=[('primary',0,30),('duplicate',1024,30),('secondary',256,30),('unmapped',4,30),('low',0,29)]
            (d/'sam').write_text(''.join(f'{q}\t{flag}\tg\t1\t{mq}\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n' for q,flag,mq in records))
            _,legacy=count_references(d/'depth',d/'sam','legacy2308')
            _,new=count_references(d/'depth',d/'sam','exclude_duplicates3332')
            self.assertEqual(legacy[0]['count'],2);self.assertEqual(new[0]['count'],1)

    def test_breadth_rounding_matches_original_awk_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);L=20001;covered=16000
            self.assertLess(covered/L,.8)
            (d/'depth').write_text(''.join(f'g\t{i}\t{int(i<=covered)}\n' for i in range(1,L+1)))
            (d/'sam').write_text('q\t0\tg\t1\t30\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
            metrics,counts=count_references(d/'depth',d/'sam','legacy2308')
            expected=subprocess.check_output(['awk','BEGIN {printf "%.2f", (16000/20001)*100}'],text=True)
            self.assertEqual(metrics[0]['breadth_percent_rounded'],float(expected))
            self.assertEqual(len(counts),1)

    def test_bam_manifest_requires_unique_bams(self):
        filters=read_tsv(ROOT/'metadata/filters.tsv')
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);rows=[]
            for r in filters:
                p=d/(r['filter_id']+'.bam');p.write_text('synthetic fixture')
                rows.append(dict(filter_id=r['filter_id'],bam=p.name))
            self.assertEqual(len(validate(filters,rows,d)),60)
            rows[1]['bam']=rows[0]['bam']
            with self.assertRaisesRegex(ValueError,'more than one filter'):validate(filters,rows,d)


if __name__=='__main__':unittest.main()
