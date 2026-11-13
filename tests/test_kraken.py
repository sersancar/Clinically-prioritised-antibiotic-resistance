"""Offline archive tests; no production database is downloaded."""
import hashlib
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from prepare_kraken import REQUIRED,prepare


class KrakenDownload(unittest.TestCase):
    def fixture(self,d,corrupt=False,omit=None,link=False):
        content={name:('fixture '+name).encode() for name in REQUIRED}
        manifest=''.join(hashlib.md5(data).hexdigest()+'  '+name+'\n' for name,data in content.items()).encode()
        (d/'checks.md5').write_bytes(manifest)
        with tarfile.open(d/'archive.tar.gz','w:gz') as archive:
            for name,data in content.items():
                if name==omit:continue
                if corrupt and name=='hash.k2d':data=b'corrupted data'
                entry=tarfile.TarInfo('release/'+name);entry.size=len(data)
                if link and name=='hash.k2d':entry.type=tarfile.SYMTYPE;entry.linkname='/etc/passwd'
                archive.addfile(entry,io.BytesIO(data))
        return content,hashlib.sha256(manifest).hexdigest()

    def run_prepare(self,d,digest):
        prepare((d/'archive.tar.gz').as_uri(),(d/'checks.md5').as_uri(),digest,d/'database',d/'taxonomy',d/'provenance.json')

    def test_archive_and_matching_taxonomy(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);content,digest=self.fixture(d);self.run_prepare(d,digest)
            self.assertEqual((d/'database/hash.k2d').read_bytes(),content['hash.k2d'])
            self.assertEqual((d/'taxonomy/nodes.dmp').read_bytes(),content['nodes.dmp'])
            self.assertFalse((d/'database/nodes.dmp').exists())
            self.assertEqual(set(json.loads((d/'provenance.json').read_text())['files']),set(REQUIRED))

    def test_corruption_is_rejected_and_outputs_removed(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);_,digest=self.fixture(d,corrupt=True)
            with self.assertRaisesRegex(ValueError,'Checksum mismatch'):self.run_prepare(d,digest)
            self.assertFalse((d/'database').exists());self.assertFalse((d/'taxonomy').exists())

    def test_changed_manifest_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);self.fixture(d)
            with self.assertRaisesRegex(ValueError,'pinned SHA256'):self.run_prepare(d,'0'*64)
            self.assertFalse((d/'database').exists())

    def test_missing_taxonomy_and_symlinks_are_rejected(self):
        for kwargs,reason in [({'omit':'nodes.dmp'},'missing required'),({'link':True},'Unsafe')]:
            with self.subTest(kwargs=kwargs),tempfile.TemporaryDirectory() as tmp:
                d=Path(tmp);_,digest=self.fixture(d,**kwargs)
                with self.assertRaisesRegex(ValueError,reason):self.run_prepare(d,digest)
                self.assertFalse((d/'database').exists())
