#!/usr/bin/env python3
import argparse
import hashlib
import json
import re
from pathlib import Path

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for x in ['lock','carpdm','crass','output']:p.add_argument('--'+x,required=True)
    a=p.parse_args();lock=json.loads(Path(a.lock).read_text());out={}
    for key in ['carpdm','crass']:
        expected=lock[key]['sha256']
        if not re.fullmatch('[a-f0-9]{64}',expected):raise ValueError(f'{key}: SHA256 must be recorded before running')
        h=hashlib.sha256()
        with open(getattr(a,key),'rb') as f:
            while b:=f.read(8*1024*1024):h.update(b)
        if h.hexdigest()!=expected:raise ValueError(f'{key}: reference checksum mismatch')
        out[key]={**lock[key],'observed_sha256':h.hexdigest(),'filename':Path(getattr(a,key)).name}
    Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
