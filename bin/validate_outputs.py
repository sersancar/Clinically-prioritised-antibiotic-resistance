#!/usr/bin/env python3
"""Fail rather than mark incomplete report generation as successful."""
import json
import sys
from pathlib import Path

if __name__=='__main__':
    root=Path(sys.argv[1]);spec=json.loads(Path(sys.argv[2]).read_text())
    missing=[x for x in spec['required'] if not (root/x).is_file() or (root/x).stat().st_size==0]
    if missing:raise SystemExit('Missing or empty manuscript outputs: '+', '.join(missing))
    # Existence is an output contract, not scientific validation. Regression
    # comparisons to the published tables remain a separate required release gate.
