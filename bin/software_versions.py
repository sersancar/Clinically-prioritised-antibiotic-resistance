#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
from pathlib import Path

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('tools',nargs='+');a=p.parse_args()
    result={'python':sys.version.split()[0]}
    for tool in a.tools:
        args=[tool] if tool=='bwa' else [tool,'-v'] if tool=='bracken' else [tool,'--version']
        r=subprocess.run(args,capture_output=True,text=True)
        lines=[line.strip() for line in (r.stdout+r.stderr).splitlines() if line.strip()]
        version=next((line for line in lines if 'Version' in line or 'version' in line),lines[0] if lines else 'unavailable')
        result[tool]=version
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
