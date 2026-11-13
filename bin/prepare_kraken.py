#!/usr/bin/env python3
"""Download a frozen Kraken2/Bracken archive and its matching taxonomy."""
import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

DATABASE=('hash.k2d','opts.k2d','taxo.k2d','database150mers.kmer_distrib')
TAXONOMY=('nodes.dmp','names.dmp')
REQUIRED=DATABASE+TAXONOMY


def parse_checksums(content):
    checks={}
    for line in content.decode().splitlines():
        if not line.strip():continue
        digest,name=line.split(maxsplit=1);name=Path(name.lstrip('*')).name
        if len(digest)!=32 or any(c not in '0123456789abcdefABCDEF' for c in digest):
            raise ValueError('Invalid published MD5 manifest')
        if name in checks:raise ValueError(f'Duplicate manifest entry: {name}')
        checks[name]=digest.lower()
    for name in REQUIRED:
        if name not in checks:raise ValueError(f'Published manifest missing {name}')
    return checks


def extract_verified(stream,checks,database,taxonomy):
    """Stream the archive, saving only regular files required by this workflow."""
    database=Path(database);taxonomy=Path(taxonomy)
    database.mkdir(parents=True,exist_ok=True);taxonomy.mkdir(parents=True,exist_ok=True)
    actual={}
    with tarfile.open(fileobj=stream,mode='r|gz') as archive:
        for member in archive:
            name=Path(member.name).name
            if name not in REQUIRED:continue
            if not member.isfile() or name in actual:raise ValueError(f'Unsafe or duplicate archive member: {member.name}')
            dest=(database if name in DATABASE else taxonomy)/name
            md5=hashlib.md5();sha=hashlib.sha256()
            with archive.extractfile(member) as src, dest.open('wb') as dst:
                while block:=src.read(8*1024*1024):
                    md5.update(block);sha.update(block);dst.write(block)
            if md5.hexdigest()!=checks[name]:
                dest.unlink();raise ValueError(f'Checksum mismatch: {name}')
            actual[name]={'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'bytes':member.size}
            print(f'Verified {name} ({member.size} bytes)',flush=True)
    missing=set(REQUIRED)-set(actual)
    if missing:raise ValueError(f'Archive missing required files: {sorted(missing)}')
    return actual


def prepare(url,checksums_url,manifest_sha256,database,taxonomy,provenance):
    with urllib.request.urlopen(checksums_url,timeout=120) as response:content=response.read()
    digest=hashlib.sha256(content).hexdigest()
    if digest!=manifest_sha256:raise ValueError('Published checksum manifest differs from the pinned SHA256')
    checks=parse_checksums(content)
    for path in [database,taxonomy]:
        if Path(path).exists():raise ValueError(f'Output already exists: {path}')
    try:
        with urllib.request.urlopen(url,timeout=120) as response:
            actual=extract_verified(response,checks,database,taxonomy)
    except Exception:
        for path in [database,taxonomy]:shutil.rmtree(path,ignore_errors=True)
        raise
    Path(provenance).write_text(json.dumps(dict(archive_url=url,checksums_url=checksums_url,
        manifest_sha256=digest,files=actual,taxonomy_source='Same archive as Kraken2/Bracken database'),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['url','checksums-url','manifest-sha256','database','taxonomy','provenance']:p.add_argument('--'+name,required=True)
    a=p.parse_args();prepare(a.url,a.checksums_url,a.manifest_sha256,a.database,a.taxonomy,a.provenance)
