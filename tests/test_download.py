"""Exercise actual curl against a local HTTP server, with no ENA downloads."""
import hashlib
import shutil
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from download_fastq import download


@unittest.skipUnless(shutil.which('curl'),'Requires curl on PATH')
class DownloadCurl(unittest.TestCase):
    def setUp(self):
        self.payload=b'ACGT verified FASTQ bytes\n'*1000
        self.requests=[];self.interrupt=False;self.corrupt=False;self.ignore_range=False
        owner=self
        class Handler(BaseHTTPRequestHandler):
            def log_message(self,*args):pass
            def do_GET(self):
                owner.requests.append((self.path,self.headers.get('Range')))
                if self.path=='/forbidden':self.send_error(403);return
                payload=b'corrupt' if owner.corrupt else owner.payload
                header=self.headers.get('Range');start=int(header.removeprefix('bytes=').split('-')[0]) if header else 0
                if owner.ignore_range:start=0
                body=payload[start:]
                self.send_response(206 if start else 200)
                self.send_header('Content-Length',str(len(body)))
                if start:self.send_header('Content-Range',f'bytes {start}-{len(payload)-1}/{len(payload)}')
                self.end_headers()
                if owner.interrupt and len(owner.requests)==1:
                    self.wfile.write(body[:100]);self.wfile.flush();self.close_connection=True
                else:
                    try:self.wfile.write(body)
                    except BrokenPipeError:pass
        self.server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
        self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
        self.url=f'http://127.0.0.1:{self.server.server_port}/read'
        self.md5=hashlib.md5(self.payload).hexdigest()
        self.tmp=tempfile.TemporaryDirectory();self.dest=Path(self.tmp.name)/'read.gz'

    def tearDown(self):
        self.server.shutdown();self.server.server_close();self.thread.join();self.tmp.cleanup()

    def test_complete_transfer_and_verified_reuse(self):
        download(self.url,self.md5,self.dest,retry_delay=0)
        self.assertEqual(self.dest.read_bytes(),self.payload)
        download(self.url,self.md5,self.dest,retry_delay=0)
        self.assertEqual(len(self.requests),1)

    def test_interrupted_transfer_resumes_at_received_offset(self):
        self.interrupt=True
        download(self.url,self.md5,self.dest,retry_delay=0)
        self.assertEqual(self.dest.read_bytes(),self.payload)
        self.assertEqual(self.requests[1][1],'bytes=100-')

    def test_existing_partial_resumes_across_invocations(self):
        self.dest.with_suffix('.gz.partial').write_bytes(self.payload[:400])
        download(self.url,self.md5,self.dest,retry_delay=0)
        self.assertEqual(self.requests[0][1],'bytes=400-')
        self.assertEqual(self.dest.read_bytes(),self.payload)

    def test_server_without_range_support_restarts_cleanly(self):
        self.ignore_range=True;self.dest.with_suffix('.gz.partial').write_bytes(self.payload[:400])
        download(self.url,self.md5,self.dest,retry_delay=0)
        self.assertEqual(self.dest.read_bytes(),self.payload)
        self.assertEqual(len(self.requests),2)

    def test_forbidden_transport_falls_back(self):
        with patch('download_fastq.transports',return_value=[self.url.replace('/read','/forbidden'),self.url]):
            download(self.url,self.md5,self.dest,retry_delay=0)
        self.assertEqual(self.dest.read_bytes(),self.payload)
        self.assertEqual(len(self.requests),2)

    def test_bad_checksum_never_replaces_existing_destination(self):
        self.corrupt=True;self.dest.write_bytes(b'previous bytes')
        with self.assertRaisesRegex(ValueError,'No verified download'):
            download(self.url,self.md5,self.dest,attempts=2,retry_delay=0)
        self.assertEqual(self.dest.read_bytes(),b'previous bytes')
        self.assertFalse(self.dest.with_suffix('.gz.partial').exists())
