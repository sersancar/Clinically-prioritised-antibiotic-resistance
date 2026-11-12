#!/bin/bash
#SBATCH --job-name=Bracken_species_relAbund_matrix_vA
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/ssanchez/AREST_scripts/logs/%x_%A.out
#SBATCH --error=/home/ssanchez/AREST_scripts/logs/%x_%A.err

set -euo pipefail

########################
# USER CONFIGURATION
########################

SAMPLE_LIST_FILE="/home/ssanchez/AREST_scripts/sample_names.txt"
BRACKEN_DIR="/data3/ssanchez/Bracken_taxonomy"
BRACKEN_LEVEL="S"

# Outputs (fractions = Option A)
MATRIX_FRAC_TSV="${BRACKEN_DIR}/Bracken_${BRACKEN_LEVEL}_fraction_matrix.tsv"
MATRIX_FULLTAX_FRAC_CSV="${BRACKEN_DIR}/Bracken_${BRACKEN_LEVEL}_fullTax_fraction_matrix.csv"

# Optional QC outputs (reads = new_est_reads)
MATRIX_READS_TSV="${BRACKEN_DIR}/Bracken_${BRACKEN_LEVEL}_newEstReads_matrix.tsv"
MATRIX_FULLTAX_READS_CSV="${BRACKEN_DIR}/Bracken_${BRACKEN_LEVEL}_fullTax_newEstReads_matrix.csv"

# TaxonKit DB (must contain nodes.dmp + names.dmp)
export TAXONKIT_DB="/data3/ssanchez/kraken2DB"

source /data3/ssanchez/miniforge3/etc/profile.d/conda.sh
# conda activate <your_env_with_taxonkit>

command -v python   >/dev/null 2>&1 || { echo "ERROR: python not in PATH" >&2; exit 1; }
command -v taxonkit >/dev/null 2>&1 || { echo "ERROR: taxonkit not in PATH" >&2; exit 1; }

while [ ! -d "/data3/ssanchez" ]; do
    echo "Directory /data3/ssanchez does not exist. Waiting..."
    sleep 5
done

echo "Job started on $(date)"
echo "Using Bracken dir    : ${BRACKEN_DIR}"
echo "Sample list file     : ${SAMPLE_LIST_FILE}"
echo "Bracken level        : ${BRACKEN_LEVEL}"
echo "Output frac TSV      : ${MATRIX_FRAC_TSV}"
echo "Output frac fullTax  : ${MATRIX_FULLTAX_FRAC_CSV}"
echo "Output reads TSV     : ${MATRIX_READS_TSV}"
echo "Output reads fullTax : ${MATRIX_FULLTAX_READS_CSV}"
echo "TAXONKIT_DB          : ${TAXONKIT_DB}"

if [[ ! -f "${SAMPLE_LIST_FILE}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

python - "${SAMPLE_LIST_FILE}" "${BRACKEN_DIR}" "${BRACKEN_LEVEL}" \
         "${MATRIX_FRAC_TSV}" "${MATRIX_FULLTAX_FRAC_CSV}" \
         "${MATRIX_READS_TSV}" "${MATRIX_FULLTAX_READS_CSV}" << 'PY'
import csv
import os
import sys
import subprocess

(sample_list_file, bracken_dir, level,
 out_frac_tsv, out_fulltax_frac_csv,
 out_reads_tsv, out_fulltax_reads_csv) = sys.argv[1:]

# ----------------------------
# 1) Collect samples & raw Bracken files
# ----------------------------
samples = []
files = []

with open(sample_list_file) as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = os.path.join(bracken_dir, f"{s}.bracken.{level}.tsv")
        if os.path.isfile(p) and os.path.getsize(p) > 0:
            samples.append(s)
            files.append(p)
        else:
            sys.stderr.write(f"WARNING: missing/empty Bracken raw file for sample '{s}': {p}\n")

ns = len(samples)
if ns == 0:
    sys.stderr.write("ERROR: No usable Bracken raw files found. Aborting.\n")
    sys.exit(1)

sys.stderr.write(f"INFO: Samples included: {ns}\n")

# ----------------------------
# 2) Build TaxID universe + matrices
#    Bracken raw format:
#      0 name
#      1 taxid
#      2 taxonomy_lvl
#      3 kraken_assigned_reads
#      4 added_reads
#      5 new_est_reads
#      6 fraction_total_reads
# ----------------------------
taxids = set()
name_by_taxid = {}
level_by_taxid = {}

# (taxid, sample_idx) -> value
mat_reads = {}
mat_frac  = {}

for idx, (s, fp) in enumerate(zip(samples, files), start=1):
    with open(fp, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r, None)
        for row in r:
            if len(row) < 7:
                continue

            taxid = row[1].strip()
            if not taxid or taxid in ("0", "1"):
                continue

            taxids.add(taxid)

            if taxid not in name_by_taxid:
                name_by_taxid[taxid] = row[0].strip()
            if taxid not in level_by_taxid:
                level_by_taxid[taxid] = row[2].strip()

            try:
                est = float(row[5])
            except ValueError:
                est = 0.0
            mat_reads[(taxid, idx)] = est

            try:
                frac = float(row[6])
            except ValueError:
                frac = 0.0
            mat_frac[(taxid, idx)] = frac

taxids_sorted = sorted(taxids, key=lambda x: int(x))

if not taxids_sorted:
    # Write empty outputs with headers
    def write_empty_tsv(path):
        with open(path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["taxid", *samples])

    def write_empty_csv(path):
        with open(path, "w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["taxid","name","level","Domain","phylum","class","order","family","genus","species", *samples])

    write_empty_tsv(out_frac_tsv)
    write_empty_tsv(out_reads_tsv)
    write_empty_csv(out_fulltax_frac_csv)
    write_empty_csv(out_fulltax_reads_csv)
    sys.stderr.write("WARNING: No taxa found in any Bracken file. Wrote empty outputs.\n")
    sys.exit(0)

# ----------------------------
# 2b) Write TaxID-only matrices
# ----------------------------
with open(out_frac_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["taxid", *samples])
    for tid in taxids_sorted:
        row = [tid] + [mat_frac.get((tid, i), 0.0) for i in range(1, ns+1)]
        w.writerow(row)

with open(out_reads_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["taxid", *samples])
    for tid in taxids_sorted:
        vals = []
        for i in range(1, ns+1):
            v = mat_reads.get((tid, i), 0.0)
            if abs(v - round(v)) < 1e-9:
                vals.append(int(round(v)))
            else:
                vals.append(v)
        w.writerow([tid, *vals])

sys.stderr.write(f"INFO: Wrote fraction matrix TSV: {out_frac_tsv}\n")
sys.stderr.write(f"INFO: Wrote new_est_reads matrix TSV: {out_reads_tsv}\n")

# ----------------------------
# 3) Add canonical ranks using TaxonKit reformat (NOT lineage-by-position)
#
# Findings applied:
# - Use {d} (Domain/superkingdom) instead of {k} (Kingdom), because {k} is often blank in NCBI.
# - Do NOT use positional parsing of lineage strings.
# - Ensure taxid is correctly associated to rank fields.
# ----------------------------
cmd = ["taxonkit", "reformat", "-I", "1", "-f", "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}"]
inp = "\n".join(taxids_sorted) + "\n"
proc = subprocess.run(cmd, input=inp, text=True, capture_output=True)

if proc.returncode != 0:
    sys.stderr.write("ERROR: taxonkit reformat failed.\n")
    sys.stderr.write(proc.stderr + "\n")
    sys.exit(proc.returncode)

taxmap = {}  # taxid -> [Domain, phylum, class, order, family, genus, species]

# TaxonKit versions differ: some print "taxid + formatted fields", others may print only formatted fields.
# We handle both robustly:
out_lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
if not out_lines:
    sys.stderr.write("WARNING: taxonkit reformat produced no output; taxonomy columns will be empty.\n")
else:
    # Case A: 8 columns => taxid + 7 ranks
    # Case B: 7 columns => only 7 ranks, assume order matches input taxids_sorted
    first = out_lines[0].split("\t")
    if len(first) >= 8:
        for line in out_lines:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 8:
                tid = parts[0]
                taxmap[tid] = parts[1:8]
    elif len(first) == 7:
        if len(out_lines) != len(taxids_sorted):
            sys.stderr.write("WARNING: taxonkit output length mismatch; taxonomy columns may be misassigned.\n")
        for tid, line in zip(taxids_sorted, out_lines):
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 7:
                taxmap[tid] = parts
    else:
        sys.stderr.write("WARNING: unexpected taxonkit output format; taxonomy columns will be empty.\n")

# ----------------------------
# 4) Write final fullTax CSVs
# ----------------------------
def write_fulltax_csv(path, mat, value_kind):
    with open(path, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "taxid", "name", "level",
            "Domain", "phylum", "class", "order", "family", "genus", "species",
            *samples
        ])
        for tid in taxids_sorted:
            domain, p, c, o, f_, g, s_ = taxmap.get(tid, ["", "", "", "", "", "", ""])
            base = [tid, name_by_taxid.get(tid, ""), level_by_taxid.get(tid, ""),
                    domain, p, c, o, f_, g, s_]
            vals = []
            for i in range(1, ns+1):
                v = mat.get((tid, i), 0.0)
                if value_kind == "reads":
                    if abs(v - round(v)) < 1e-9:
                        vals.append(int(round(v)))
                    else:
                        vals.append(v)
                else:
                    vals.append(v)
            w.writerow(base + vals)

write_fulltax_csv(out_fulltax_frac_csv, mat_frac, "frac")
write_fulltax_csv(out_fulltax_reads_csv, mat_reads, "reads")

sys.stderr.write(f"INFO: Wrote full taxonomy fraction CSV: {out_fulltax_frac_csv}\n")
sys.stderr.write(f"INFO: Wrote full taxonomy reads CSV: {out_fulltax_reads_csv}\n")
PY

echo "[$(date)] Finished."
echo "Fractions (TaxID TSV)       : ${MATRIX_FRAC_TSV}"
echo "Fractions (FullTax CSV)     : ${MATRIX_FULLTAX_FRAC_CSV}"
echo "Reads new_est (TaxID TSV)   : ${MATRIX_READS_TSV}"
echo "Reads new_est (FullTax CSV) : ${MATRIX_FULLTAX_READS_CSV}"

