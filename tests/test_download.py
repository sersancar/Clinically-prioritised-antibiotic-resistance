import hashlib
import io
import ssl
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

    def test_ssl_failure_mid_stream_falls_back_to_http(self):
        class BrokenTLS(io.BytesIO):
            def read(self,n=-1):
                if self.tell():raise ssl.SSLError('[SYS] unknown error (_ssl.c:2559)')
                return super().read(4)
        payload=b'verified complete read file'
        url='https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR1/ERR1_2.fastq.gz'
        with tempfile.TemporaryDirectory() as tmp,patch('download_fastq.urllib.request.urlopen',side_effect=[BrokenTLS(b'partial bytes'),io.BytesIO(payload)]) as opener:
            dest=Path(tmp)/'read.gz'
            download(url,hashlib.md5(payload).hexdigest(),dest)
            self.assertEqual(dest.read_bytes(),payload)
            self.assertEqual(opener.call_count,2)
            self.assertTrue(opener.call_args_list[1].args[0].full_url.startswith('http://'))
            self.assertFalse(dest.with_suffix('.gz.partial').exists())

    def test_checksum_failure_falls_back_and_verifies_replacement(self):
        payload=b'verified complete read file'
        url='https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR1/ERR1_2.fastq.gz'
        with tempfile.TemporaryDirectory() as tmp,patch('download_fastq.urllib.request.urlopen',side_effect=[io.BytesIO(b'truncated'),io.BytesIO(payload)]) as opener:
            dest=Path(tmp)/'read.gz'
            download(url,hashlib.md5(payload).hexdigest(),dest)
            self.assertEqual(dest.read_bytes(),payload)
            self.assertTrue(opener.call_args_list[1].args[0].full_url.startswith('http://'))
