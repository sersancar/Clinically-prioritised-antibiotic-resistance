import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from ena_manifest import read_tsv
from prepare_r_inputs import legacy_id,prepare


class OriginalIntegration(unittest.TestCase):
    def test_archived_sources_unchanged(self):
        for record in json.loads((ROOT/'provenance/original_scripts/manifest.json').read_text()):
            self.assertEqual(hashlib.sha256((ROOT/record['archived_path']).read_bytes()).hexdigest(),record['sha256'])

    def test_filter_matches_original_scoring(self):
        source=(ROOT/'provenance/original_scripts/kraken2Filter.sh').read_text()
        start=source.index("<< 'EOF'\n")+len("<< 'EOF'\n")
        original=source[start:source.index('\nEOF',start)]
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'original.py').write_text(original)
            raw=('C\tboundary\tA taxon (taxid 10)\t150|150\t10:2 20:8 0:99 A:80\n'
                 'C\tlow\tA taxon (taxid 10)\t150|150\t10:1 20:9\n'
                 'C\tunknown\t10\t150|150\t0:10 A:5\n'
                 'C\tempty\t10\t150|150\t\n'
                 'C\tboundary\tA taxon (taxid 10)\t150|150\t10:2 20:8\n'
                 'U\tno\t0\t150|150\t0:20\n')
            (d/'input').write_text(raw)
            for executable,prefix in [(d/'original.py','original'),(ROOT/'bin/kraken_score_filter.py','portable')]:
                subprocess.run([sys.executable,str(executable),'.2',str(d/'input'),str(d/(prefix+'.out')),str(d/(prefix+'.map'))],check=True)
            self.assertEqual((d/'original.out').read_text(),(d/'portable.out').read_text())
            self.assertEqual(set((d/'original.map').read_text().splitlines()),set((d/'portable.map').read_text().splitlines()))
            self.assertNotIn('\tlow\t',(d/'portable.out').read_text())

    def test_report_counts_clades_and_preserves_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/'nodes.dmp').write_text('1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n10\t|\t2\t|\tgenus\t|\n11\t|\t10\t|\tspecies\t|\n')
            (d/'names.dmp').write_text(''.join(f'{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n' for tid,name in [(1,'root'),(2,'Bacteria'),(10,'Test genus'),(11,'Test species')]))
            (d/'filtered').write_text('C\ta\tTest genus (taxid 10)\t150|150\t10:2\nC\tb\tTest species (taxid 11)\t150|150\t11:2\n')
            subprocess.run([sys.executable,str(ROOT/'bin/build_kraken_report.py'),str(d/'filtered'),str(d),str(d/'report')],check=True)
            rows=[r.split('\t') for r in (d/'report').read_text().splitlines()]
            byid={r[4]:r for r in rows}
            self.assertEqual(byid['1'][1],'2');self.assertEqual(byid['10'][1],'2')
            self.assertEqual(byid['10'][2],'1');self.assertEqual(byid['11'][5].strip(),'Test species')

    def test_portable_r_metadata_uses_runtime_totals(self):
        filters=read_tsv(ROOT/'metadata/filters.tsv')
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            for name in ['amr','metrics','tax','qc']:(d/name).mkdir()
            ref='MEG_1|Drugs|MLS|Mechanism|gene'
            (d/'megares.fa').write_text('>'+ref+'\nAAAA\n')
            for r in filters:
                fid=r['filter_id'];(d/'amr'/(fid+'.amr')).mkdir()
                (d/'amr'/(fid+'.amr')/'reference_counts.tsv').write_text('gene_accession\tcount\n'+ref+'\t2\n')
                (d/'metrics'/(fid+'.carpdm.metrics.tsv')).write_text('filter_id\treference\tlength\tcovered_bases\tbreadth\tmean_depth\tmapped_read_ends\n'+fid+'\tg\t4\t4\t1\t5\t2\n')
                (d/'tax'/(fid+'.CARD_ARG_taxonomy.tsv')).write_text('ARG_ID\tread_count\tgenus\n')
                (d/'qc'/(fid+'.fastp.json')).write_text(json.dumps({'summary':{'after_filtering':{'total_reads':100}}}))
            prepare(filters,d/'amr',d/'metrics',d/'tax',d/'qc',d/'megares.fa',d/'out')
            totals=read_tsv(d/'out/total_reads.tsv')
            self.assertEqual(len(totals),60);self.assertEqual({r['total_reads'] for r in totals},{'100'})
            self.assertEqual(legacy_id(next(r for r in filters if r['filter_id']=='D_1_290321')),'F1FW290321')
            self.assertEqual(len(read_tsv(d/'out/sample_manifest.tsv')),60)


if __name__=='__main__':unittest.main()
