#!/usr/bin/env python3
"""Stream an ENA FASTQ to disk, atomically retaining only verified bytes."""
import argparse
import hashlib
import os
import urllib.request
from pathlib import Path


def download(url, expected, destination):
    destination = Path(destination)
    temporary = destination.with_suffix(destination.suffix+'.partial')
    digest = hashlib.md5()
    try:
        with urllib.request.urlopen(url, timeout=120) as response, temporary.open('wb') as out:
            while chunk := response.read(8*1024*1024):
                digest.update(chunk);out.write(chunk)
        if digest.hexdigest().lower() != expected.lower():
            raise ValueError(f'MD5 mismatch for {destination.name}')
        os.replace(temporary,destination)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('url');p.add_argument('md5');p.add_argument('destination')
    a=p.parse_args();download(a.url,a.md5,a.destination)
