import csv
import hashlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from ena_manifest import read_tsv,resolve,validate_filters
from coverage_metrics import summarise
from pool_screening import pool
from download_fastq import download
from map_ena_titles import canonical,map_titles


class Methods(unittest.TestCase):
    def setUp(self):
        self.filters=read_tsv(ROOT/'metadata/filters.tsv')

    def test_study_design(self):
        groups=validate_filters(self.filters,require_runs=False)
        self.assertEqual(len(groups),20)
        self.assertEqual(sum(r['site']=='C' for r in self.filters),6)
        for group in groups.values():
            self.assertEqual(len(group),3)

    def test_mapping_is_required(self):
        rows=[dict(r,run_accession='') for r in self.filters]
        with self.assertRaisesRegex(ValueError,'mapping missing'):
            validate_filters(rows)

    def test_mapping_titles_bijection(self):
        rows=map_titles(self.filters,read_tsv(ROOT/'metadata/ena_PRJEB75359_2026-11-12.tsv'))
        self.assertEqual([r['run_accession'] for r in rows],[r['run_accession'] for r in self.filters])
        self.assertEqual(canonical('1Freshwater140421'),'D_1_140421')
        with self.assertRaisesRegex(ValueError,'Unrecognised'):canonical('34')

    def test_run_cannot_belong_to_two_filters(self):
        rows=[dict(r,run_accession=f'ERR{i+1}') for i,r in enumerate(self.filters)]
        rows[1]['run_accession']=rows[0]['run_accession']
        with self.assertRaisesRegex(ValueError,'multiple filters'):
            validate_filters(rows)

    def test_ena_checks_paired_order_and_study_membership(self):
        rows=[dict(r,run_accession=f'ERR{i+1}') for i,r in enumerate(self.filters)]
        ena=[dict(run_accession=r['run_accession'],library_layout='PAIRED',
                  fastq_ftp=f"ftp.sra.ebi.ac.uk/{r['run_accession']}_1.fastq.gz;ftp.sra.ebi.ac.uk/{r['run_accession']}_2.fastq.gz",
                  fastq_md5='0'*32+';'+'1'*32) for r in rows]
        self.assertEqual(len(resolve(rows,ena)),60)
        ena[0]['library_layout']='SINGLE'
        with self.assertRaisesRegex(ValueError,'exactly two'):resolve(rows,ena)

    def test_full_length_denominator_and_read_end_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'ref').write_text('>gene\nAAAAAAAAAA\n>zero\nTTTT\n')
            # Five overlapping read ends each cover the same eight bases.
            (d/'depth').write_text(''.join(f'gene\t{i}\t5\n' for i in range(1,9)))
            (d/'sam').write_text(''.join(f'q{i}\t0\tgene\t1\t30\t8M\t*\t0\t0\tAAAAAAAA\tIIIIIIII\n' for i in range(5)))
            rows=summarise(d/'ref',d/'depth',d/'sam')
            self.assertEqual(rows[0]['breadth'],.8)
            self.assertEqual(rows[0]['mean_depth'],4)
            self.assertEqual(rows[0]['mapped_read_ends'],5)
            self.assertFalse(rows[0]['primary']);self.assertTrue(rows[0]['relaxed'])
            self.assertEqual(rows[1]['breadth'],0)

    def test_pool_keeps_failed_filters_in_denominator(self):
        totals={r['filter_id']:100 for r in self.filters}
        metrics=[dict(filter_id=r['filter_id'],reference='g',length=1000,
                      mapped_read_ends=30 if r['replicate']=='1' else 200,
                      primary=r['replicate']=='1',relaxed=True) for r in self.filters]
        rows=pool(self.filters,metrics,totals)
        primary=next(r for r in rows if r['screen']=='primary')
        self.assertEqual(primary['counts'],30)
        self.assertEqual(primary['total_read_ends'],300)
        self.assertEqual(primary['RPM'],100000)
        self.assertEqual(primary['RPKM'],100000)
        self.assertEqual(primary['positive_filters'],1)

    def test_missing_filter_metrics_fail(self):
        totals={r['filter_id']:100 for r in self.filters}
        metrics=[dict(filter_id=r['filter_id'],reference='g',length=1000,
                      mapped_read_ends=0,primary=False,relaxed=False) for r in self.filters[:-1]]
        with self.assertRaisesRegex(ValueError,'Incomplete'):pool(self.filters,metrics,totals)

    def test_checksum_failure_removes_partial_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'r.fastq.gz'
            with patch('urllib.request.urlopen',return_value=io.BytesIO(b'corrupt')):
                with self.assertRaisesRegex(ValueError,'MD5 mismatch'):download('unused','0'*32,dest)
            self.assertFalse(dest.exists());self.assertFalse(dest.with_suffix('.gz.partial').exists())

    def test_verified_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'r.fastq.gz';data=b'compressed-bytes'
            with patch('urllib.request.urlopen',return_value=io.BytesIO(data)):
                download('unused',hashlib.md5(data).hexdigest(),dest)
            self.assertEqual(dest.read_bytes(),data)


if __name__=='__main__':unittest.main()
