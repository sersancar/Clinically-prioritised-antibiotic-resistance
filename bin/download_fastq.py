#!/usr/bin/env python3
"""Download ENA FASTQs with transport fallback and mandatory MD5 verification."""
import argparse
import hashlib
import http.client
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def checksum(path):
    digest=hashlib.md5()
    with Path(path).open('rb') as stream:
        while chunk:=stream.read(8*1024*1024):digest.update(chunk)
    return digest.hexdigest()


def transports(url):
    parsed=urllib.parse.urlsplit(url)
    if parsed.hostname=='ftp.sra.ebi.ac.uk' and parsed.scheme in ['https','http','ftp']:
        schemes=[parsed.scheme]+[s for s in ['https','http','ftp'] if s!=parsed.scheme]
        return [urllib.parse.urlunsplit(parsed._replace(scheme=s)) for s in schemes]
    return [url]


def download(url, expected, destination):
    destination=Path(destination)
    if destination.is_file() and checksum(destination).lower()==expected.lower():
        print(f'Reusing MD5-verified {destination.name}',flush=True)
        return
    temporary=destination.with_suffix(destination.suffix+'.partial')
    failures=[]
    for address in transports(url):
        for attempt in range(2):
            try:
                # A Request sets a descriptive User-Agent, but changes neither file
                # bytes nor the ENA URL/checksum crosswalk produced by the resolver.
                request=urllib.request.Request(address,headers={'User-Agent':'coastal-amr-nextflow/0.2'}) if urllib.parse.urlsplit(address).scheme else address
                digest=hashlib.md5()
                print(f'Downloading {destination.name} using {address}',flush=True)
                with urllib.request.urlopen(request,timeout=120) as response,temporary.open('wb') as out:
                    while chunk:=response.read(8*1024*1024):
                        digest.update(chunk);out.write(chunk)
                if digest.hexdigest().lower()!=expected.lower():
                    raise ValueError(f'MD5 mismatch for {destination.name}: expected {expected.lower()}, got {digest.hexdigest()}')
                os.replace(temporary,destination)
                print(f'Verified {destination.name}',flush=True)
                return
            except (urllib.error.URLError,TimeoutError,ConnectionError,ssl.SSLError,http.client.HTTPException,ValueError) as exc:
                failures.append(f'{address}: {type(exc).__name__}: {exc}')
                print(f'Download attempt failed: {failures[-1]}',file=sys.stderr,flush=True)
                # Protocol/permission errors need another transport, not an
                # immediate repeat of the same forbidden URL.
                if isinstance(exc,(ssl.SSLError,ValueError)):break
                if isinstance(exc,urllib.error.HTTPError) and exc.code not in [408,429,500,502,503,504]:break
                if attempt==0:time.sleep(2)
            finally:
                temporary.unlink(missing_ok=True)
    raise ValueError('No verified download could be obtained:\n'+'\n'.join(failures))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('url');p.add_argument('md5');p.add_argument('destination')
    a=p.parse_args();download(a.url,a.md5,a.destination)
