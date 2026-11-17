import hashlib
import io
import sys
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from download_fastq import download


class DownloadFallback(unittest.TestCase):
    def test_forbidden_https_falls_back_to_http(self):
        payload=b'valid compressed file bytes';md5=hashlib.md5(payload).hexdigest()
        url='https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR1/ERR1_1.fastq.gz'
        forbidden=urllib.error.HTTPError(url,403,'Forbidden',{},None)
        with tempfile.TemporaryDirectory() as tmp,patch('download_fastq.urllib.request.urlopen',side_effect=[forbidden,io.BytesIO(payload)]) as opener:
            dest=Path(tmp)/'read.gz';download(url,md5,dest)
            self.assertEqual(dest.read_bytes(),payload)
            self.assertEqual(opener.call_args_list[1].args[0].full_url,url.replace('https://','http://'))
            self.assertFalse(dest.with_suffix('.gz.partial').exists())

    def test_verified_mate_is_not_downloaded_again(self):
        with tempfile.TemporaryDirectory() as tmp,patch('download_fastq.urllib.request.urlopen') as opener:
            dest=Path(tmp)/'read.gz';dest.write_bytes(b'verified')
            download('https://ftp.sra.ebi.ac.uk/x',hashlib.md5(b'verified').hexdigest(),dest)
            opener.assert_not_called()

    def test_corruption_never_replaces_an_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp,patch('download_fastq.urllib.request.urlopen',side_effect=lambda *a,**k:io.BytesIO(b'corrupt')),patch('download_fastq.time.sleep'):
            dest=Path(tmp)/'read.gz';dest.write_bytes(b'previous file')
            with self.assertRaisesRegex(ValueError,'No verified download'):download('https://example.test/read.gz','0'*32,dest)
            self.assertEqual(dest.read_bytes(),b'previous file')
            self.assertFalse(dest.with_suffix('.gz.partial').exists())
