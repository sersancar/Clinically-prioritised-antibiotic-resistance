import csv
import gzip
import io
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'));sys.path.insert(0,str(ROOT/'tests'))
from common import Taxonomy,query_id,read_tsv,write_tsv
from ena_manifest import resolve,validate_filters
from download_fastq import download
from aggregate_results import assemble
from link_reads import support_score,taxid_from_label,kraken_report
from make_fixture import make
from resources import copy_fasta


class Methods(unittest.TestCase):
    def test_reference_decompression_after_staging(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);payload=b'>ref|sul2\nACGTACGT\n'
            staged=d/'input_carpdm.fa'
            staged.write_bytes(gzip.compress(payload))
            copy_fasta(staged,d/'carpdm.fa')
            self.assertEqual((d/'carpdm.fa').read_bytes(),payload)
            staged.write_bytes(payload)
            copy_fasta(staged,d/'plain.fa')
            self.assertEqual((d/'plain.fa').read_bytes(),payload)

    def test_design_and_mapping(self):
        filters=read_tsv(ROOT/'data/filters.tsv');groups=validate_filters(filters)
        self.assertEqual(len(groups),20);self.assertEqual(sum(r['site']=='C' for r in filters),6)
        resolved=resolve(filters,read_tsv(ROOT/'data/ena_PRJEB75359_2026-11-12.tsv'))
        self.assertEqual(len(resolved),60)
        with self.assertRaises(ValueError):validate_filters(filters[:-1])

    def test_checksum_atomic_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'reads.gz'
            with patch('urllib.request.urlopen',return_value=io.BytesIO(b'corrupt')):
                with self.assertRaises(ValueError):download('unused','0'*32,dest)
            self.assertFalse(dest.exists());self.assertFalse(dest.with_suffix('.gz.partial').exists())

    def test_custom_score_and_pair_identifiers(self):
        self.assertEqual(support_score(10,'10:2 20:8 A:100 0:999 |:|'),.2)
        self.assertIsNone(support_score(10,'A:100 0:999'))
        self.assertLess(support_score(10,'10:1 20:9'),.2)
        self.assertEqual(taxid_from_label('Taxon name with spaces (taxid 10)'),10)
        self.assertEqual(query_id('pair/1 comment'),'pair');self.assertEqual(query_id('pair/2'),'pair')

    def test_taxonomy_ranks_and_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/'nodes.dmp').write_text('1 | 1 | no rank |\n2 | 1 | superkingdom |\n10 | 2 | class |\n20 | 10 | genus |\n30 | 20 | species |\n')
            (d/'names.dmp').write_text(''.join(f'{i} | {name} | | scientific name |\n' for i,name in [(1,'root'),(2,'Bacteria'),(10,'A class'),(20,'A genus'),(30,'A species')]))
            tax=Taxonomy(d);info=tax.annotate(30)
            self.assertEqual(info['genus'],'A genus');self.assertEqual(info['class'],'A class');self.assertEqual(info['domain'],'Bacteria')
            kraken_report({30:2,20:1},tax,d/'report')
            rows=[r.split('\t') for r in (d/'report').read_text().splitlines()]
            self.assertEqual(next(r for r in rows if r[4]=='1')[1],'3')
            with self.assertRaises(ValueError):tax.annotate(999)

    def test_pooling_and_sensitivity(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);make(ROOT,d/'inputs')
            assemble(ROOT/'data/filters.tsv',d/'inputs',ROOT/'data',d/'data')
            self.assertFalse(list((d/'data/tables').glob('TableS[23]_*.tsv')))
            self.assertTrue((d/'data/provenance/software_versions.tsv').exists())
            rows=read_tsv(d/'data/tables/TableS5_CARPDM_biological.tsv')
            bio=[r for r in rows if r['biosample']=='A_2021-03-29']
            primary=[r for r in bio if r['screen']=='primary'];relaxed=[r for r in bio if r['screen']=='relaxed']
            self.assertTrue(all(float(r['RPM'])==0 for r in primary))
            self.assertTrue(all(int(r['positive_filters'])==1 for r in relaxed))
            expected_den=sum(1000000+rep*100 for rep in [1,2,3])
            self.assertEqual(int(relaxed[0]['total_read_ends']),expected_den)
            self.assertAlmostEqual(float(relaxed[0]['RPM']),4e6/expected_den)
            self.assertAlmostEqual(float(relaxed[0]['RPKM']),float(relaxed[0]['RPM']))
            comp=read_tsv(d/'data/tables/TableS7_completeness.tsv')
            nohit=next(r for r in comp if r['catalogue']=='carpdm' and r['level']=='site' and r['identifier']=='A')
            self.assertEqual(nohit['genus_percent'],'')

    def test_missing_filter_is_not_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);make(ROOT,d/'inputs');(d/'inputs/A_1_290321.carpdm.metrics.tsv').unlink()
            with self.assertRaisesRegex(ValueError,'every filter'):
                assemble(ROOT/'data/filters.tsv',d/'inputs',ROOT/'data',d/'data')


try:import pysam
except ImportError:pysam=None


@unittest.skipUnless(pysam and shutil.which('samtools'),'Requires the bioinformatics environment')
class BAMMethods(unittest.TestCase):
    def test_overlapping_ends_and_full_length_depth(self):
        from bam_metrics import metrics
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'ref.fa').write_text('>g\n'+'A'*100+'\n')
            bam=d/'test.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'g','LN':100}]}) as f:
                for i in range(10):
                    r=pysam.AlignedSegment();r.query_name='pair'+str(i//2);r.query_sequence='A'*50;r.flag=99 if i%2==0 else 147
                    r.reference_id=0;r.reference_start=0;r.mapping_quality=30;r.cigar=[(0,50)]
                    r.next_reference_id=0;r.next_reference_start=0;r.template_length=50 if i%2==0 else -50
                    r.query_qualities=pysam.qualitystring_to_array('I'*50);f.write(r)
            pysam.index(str(bam));row=metrics(str(bam),d/'ref.fa','filter','carpdm')[0]
            self.assertEqual(row['mapped_read_ends'],10);self.assertEqual(row['breadth'],.5)
            self.assertEqual(row['mean_depth'],5);self.assertFalse(row['primary'])
            marker=metrics(str(bam),d/'ref.fa','filter','crass')[0]
            self.assertEqual(marker['mapped_read_ends'],10);self.assertEqual(marker['breadth'],.5)

    def test_unverifiable_snp_target_is_zero(self):
        from verify_snps import verify
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);bundle=d/'bundle';(bundle/'data').mkdir(parents=True)
            (bundle/'data/SNPinfo.fasta').write_text('>MEG_2|Drugs|MLS|M|g|Nuc:A1T\nAAAA\n')
            ref='MEG_1|Drugs|MLS|M|g|RequiresSNPConfirmation'
            write_tsv(d/'metrics.tsv',[dict(reference=ref,mapped_read_ends=10,primary=True,length=100,breadth=1,mean_depth=10)])
            verify('unused',d/'metrics.tsv',bundle,'unused.fa','sample',1,d/'output')
            self.assertEqual(read_tsv(d/'output/counts.tsv')[0]['count'],'0')

    @unittest.skipUnless(os.environ.get('AMR_SNP_BUNDLE'),'Set AMR_SNP_BUNDLE to the pinned verifier for the allele test')
    def test_known_mutant_and_wildtype_counts(self):
        from Bio import SeqIO
        from verify_snps import verify
        bundle=Path(os.environ['AMR_SNP_BUNDLE'])
        with (bundle/'data/SNPinfo.fasta').open() as handle:record=next(SeqIO.parse(handle,'fasta'))
        ref='|'.join(record.id.split('|')[:5])+'|RequiresSNPConfirmation';sequence=str(record.seq)
        variant=re.search(r'Nuc:([ACGT])(\d+)([ACGT])',record.id);base,pos,alt=variant.groups();pos=int(pos)-1
        self.assertEqual(sequence[pos],base)
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'ref.fa').write_text('>'+ref+'\n'+sequence+'\n')
            with pysam.AlignmentFile(d/'sample.bam','wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':ref,'LN':len(sequence)}]}) as bam:
                for i in range(10):
                    start=pos-60;seq=list(sequence[start:start+150]);seq[60]=alt if i<5 else base
                    r=pysam.AlignedSegment();r.query_name='read'+str(i);r.query_sequence=''.join(seq);r.flag=0;r.reference_id=0
                    r.reference_start=start;r.mapping_quality=60;r.cigar=[(0,150)];r.query_qualities=pysam.qualitystring_to_array('I'*150);bam.write(r)
            write_tsv(d/'metrics.tsv',[dict(reference=ref,mapped_read_ends=10,primary=True,length=len(sequence),breadth=1,mean_depth=8)])
            verify(d/'sample.bam',d/'metrics.tsv',bundle,d/'ref.fa','sample',1,d/'verified')
            self.assertEqual(read_tsv(d/'verified/counts.tsv')[0]['count'],'5')


if __name__=='__main__':unittest.main()
