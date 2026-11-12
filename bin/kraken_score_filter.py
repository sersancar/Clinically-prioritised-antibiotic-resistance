#!/usr/bin/env python3
# Scoring logic extracted from the Python block from provenance/original_scripts/kraken2Filter.sh
import sys, re

min_score = float(sys.argv[1])
infile = sys.argv[2]
filtered_out = sys.argv[3]
read2taxid_out = sys.argv[4]

taxid_pattern = re.compile(r"\(taxid\s+(\d+)\)")

def extract_label_taxid(name_field: str) -> str:
    s = name_field.strip()
    if s.isdigit():
        return s
    m = taxid_pattern.search(s)
    if m:
        return m.group(1)
    parts = s.split()
    return parts[-1] if parts else "0"

def compute_score(label_taxid: str, kmer_field: str):
    if not kmer_field.strip():
        return None
    C = 0
    Q = 0
    for token in kmer_field.split():
        if ":" not in token:
            continue
        taxid, count = token.split(":", 1)
        if taxid in ("A", "0"):
            continue
        try:
            n = int(count)
        except ValueError:
            continue
        Q += n
        if taxid == label_taxid:
            C += n
    if Q == 0:
        return None
    return C / Q


with open(infile) as fin, open(filtered_out, "w") as fout, open(read2taxid_out, "w") as fr2t:
    for line in fin:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        if parts[0] != "C":
            # ignore unclassified reads
            continue

        read_id = parts[1]
        name_field = parts[2]
        label_taxid = extract_label_taxid(name_field)
        kmer_field = parts[4] if len(parts) > 4 else ""

        score = compute_score(label_taxid, kmer_field)
        keep = (score is None) or (score >= min_score)
        if not keep:
            continue

        fout.write(line)

        # Deduplication is performed by external sort in the adapter.
        fr2t.write(f"{read_id}\t{label_taxid}\n")
