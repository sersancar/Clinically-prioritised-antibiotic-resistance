#!/usr/bin/env python3
"""Download ENA FASTQs with resumable curl transfers and mandatory MD5 checks."""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path


def checksum(path):
    digest=hashlib.md5()
    with Path(path).open('rb') as stream:
        while chunk:=stream.read(8*1024*1024):digest.update(chunk)
    return digest.hexdigest()


def transports(url):
    parsed=urllib.parse.urlsplit(url)
    if parsed.hostname=='ftp.sra.ebi.ac.uk' and parsed.scheme in ['https','http','ftp']:
        # HTTP was confirmed to work on the affected workstation. Every transport
        # must produce exactly the same ENA checksum before a file is accepted.
        return [urllib.parse.urlunsplit(parsed._replace(scheme=s)) for s in ['http','https','ftp']]
    return [url]


def download(url,expected,destination,attempts=3,retry_delay=2):
    expected=expected.lower()
    if len(expected)!=32 or any(c not in '0123456789abcdef' for c in expected):
        raise ValueError('Expected a valid ENA MD5 checksum')
    destination=Path(destination)
    if destination.is_file() and checksum(destination)==expected:
        print(f'Reusing MD5-verified {destination.name}',flush=True);return
    curl=shutil.which('curl')
    if not curl:raise RuntimeError('curl is required on the Nextflow launch host; install curl and resume the workflow')
    temporary=destination.with_suffix(destination.suffix+'.partial')
    # A previous invocation may have finished receiving bytes before being killed.
    if temporary.is_file() and checksum(temporary)==expected:
        os.replace(temporary,destination);print(f'Verified retained partial {destination.name}',flush=True);return
    failures=[]
    for address in transports(url):
        for attempt in range(attempts):
            offset=temporary.stat().st_size if temporary.exists() else 0
            print(f'Downloading {destination.name}: {address}; attempt {attempt+1}/{attempts}; resume offset {offset}',flush=True)
            command=[curl,'--disable','--fail','--location','--max-redirs','5',
                '--silent','--show-error','--connect-timeout','30','--speed-limit','1024',
                '--speed-time','120','--http1.1','--continue-at','-',
                '--output',str(temporary),'--url',address]
            result=subprocess.run(command,check=False)
            if result.returncode==0:
                actual=checksum(temporary) if temporary.exists() else hashlib.md5(b'').hexdigest()
                if actual==expected:
                    os.replace(temporary,destination);print(f'Verified {destination.name}',flush=True);return
                failures.append(f'{address}: MD5 mismatch: expected {expected}, got {actual}')
                # Never resume a completed transfer whose contents failed MD5.
                temporary.unlink(missing_ok=True)
            else:
                failures.append(f'{address}: curl exit status {result.returncode}; partial retained for resumption')
                if result.returncode in [33,36]:
                    # Server refused resumption: allow the next attempt to start fresh.
                    temporary.unlink(missing_ok=True)
                # Access/not-found errors need another transport, not repeated calls.
                if result.returncode in [22,67,78]:
                    print(failures[-1],file=sys.stderr,flush=True);break
                # Filesystem errors should be fixed locally rather than retried.
                if result.returncode in [23,26]:raise RuntimeError(failures[-1])
            print(failures[-1],file=sys.stderr,flush=True)
            if attempt+1<attempts:time.sleep(retry_delay*(attempt+1))
    raise ValueError('No verified download could be obtained:\n'+'\n'.join(failures))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('url');p.add_argument('md5');p.add_argument('destination')
    a=p.parse_args();download(a.url,a.md5,a.destination)
