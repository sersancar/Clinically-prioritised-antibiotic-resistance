#!/bin/bash
#SBATCH --job-name=Kraken2Filter
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-60
#SBATCH --output=/home/ssanchez/AREST_scripts/logs/%x_%A_%a.out
#SBATCH --error=/home/ssanchez/AREST_scripts/logs/%x_%A_%a.err

set -euo pipefail

########################
# USER CONFIGURATION
########################

# File with sample names (one per line, no header), 60 lines
SAMPLE_LIST_FILE="/home/ssanchez/AREST_scripts/sample_names.txt"

# Kraken2 results directory
KRAKEN_RES_DIR="/data3/ssanchez/kraken2"
KRAKEN_SUFFIX=".out"

# Output directory
OUTDIR="/data3/ssanchez/Kraken_Taxonomy"
mkdir -p "${OUTDIR}"

# Control strength of Kraken2 taxonomic assignments (score on k-mer field)
KRAKEN_SCORE=0.2

# TaxonKit DB: MUST be the directory containing nodes.dmp and names.dmp.
export TAXONKIT_DB="/data3/ssanchez/kraken2DB"

# Conda env (make sure python and taxonkit are available)
source /data3/ssanchez/miniforge3/etc/profile.d/conda.sh
# conda activate <your_env>

command -v python   >/dev/null 2>&1 || { echo "ERROR: python not in PATH" >&2; exit 1; }
command -v taxonkit >/dev/null 2>&1 || { echo "ERROR: taxonkit not in PATH" >&2; exit 1; }

# Wait until main directory is mounted / available
while [ ! -d "/data3/ssanchez" ]; do
    echo "Directory /data3/ssanchez does not exist. Waiting..."
    sleep 5
done

echo "Job started on $(date)"
echo "Running on node(s): ${SLURM_NODELIST:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "KRAKEN_RES_DIR: ${KRAKEN_RES_DIR}"
echo "OUTDIR       : ${OUTDIR}"
echo "KRAKEN_SCORE : ${KRAKEN_SCORE}"
echo "TAXONKIT_DB  : ${TAXONKIT_DB}"

########################
# GET SAMPLE NAME (ARRAY MODE; optional CLI override)
########################
if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "ERROR: No sample argument provided and SLURM_ARRAY_TASK_ID is not set." >&2
        exit 1
    fi
    if [[ ! -f "${SAMPLE_LIST_FILE}" ]]; then
        echo "ERROR: Sample list file not found: ${SAMPLE_LIST_FILE}" >&2
        exit 1
    fi
    SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST_FILE}" | tr -d '\r' | xargs)"
    if [[ -z "${SAMPLE}" ]]; then
        echo "ERROR: No sample name on line ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST_FILE}" >&2
        exit 1
    fi
fi

echo "Processing sample: ${SAMPLE}"

KRAKEN_ORIG="${KRAKEN_RES_DIR}/${SAMPLE}${KRAKEN_SUFFIX}"
if [[ ! -f "${KRAKEN_ORIG}" ]]; then
    echo "ERROR: Kraken file not found: ${KRAKEN_ORIG}" >&2
    exit 1
fi

PREFIX="${OUTDIR}/${SAMPLE}"
KRAKEN_FILTERED="${PREFIX}.kraken.score${KRAKEN_SCORE}.out"

# Outputs
READ2TAXID_TSV="${PREFIX}.read_to_taxid.tsv"
COUNTS_NUM_TSV="${PREFIX}.taxid_counts.numeric.tsv"
TAXIDS_TXT="${PREFIX}.taxids.txt"
TAXLINEAGE="${PREFIX}.taxid_lineage.tsv"
TAXRANKS="${PREFIX}.taxid_ranks.tsv"
TAXMAP="${PREFIX}.taxid_taxonomy.tsv"
FINAL_TABLE="${OUTDIR}/${SAMPLE}.taxonomy.tsv"

########################
# 1) FILTER KRAKEN BY SCORE, AND WRITE READ -> TAXID (NUMERIC)
########################
echo "[$(date)] Filtering Kraken2 output by score >= ${KRAKEN_SCORE} -> ${KRAKEN_FILTERED}"
echo "[$(date)] Writing read -> numeric taxid -> ${READ2TAXID_TSV}"

python - "$KRAKEN_SCORE" "${KRAKEN_ORIG}" "${KRAKEN_FILTERED}" "${READ2TAXID_TSV}" << 'EOF'
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

seen = set()

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

        key = (read_id, label_taxid)
        if key not in seen:
            fr2t.write(f"{read_id}\t{label_taxid}\n")
            seen.add(key)
EOF

LC_ALL=C sort -t $'\t' -k1,1 -u "${READ2TAXID_TSV}" -o "${READ2TAXID_TSV}"

N_READS=$(wc -l < "${READ2TAXID_TSV}" || echo 0)
echo "  [INFO] Classified reads after filtering: ${N_READS}"

if [[ "${N_READS}" -le 0 ]]; then
    echo "WARNING: No classified reads after filtering. Writing empty taxonomy table."
    echo -e "taxid\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
    exit 0
fi

########################
# 2) AGGREGATE: TAXID COUNTS (NUMERIC)
########################
echo "[$(date)] Aggregating read counts per numeric taxid -> ${COUNTS_NUM_TSV}"

awk -F'\t' 'BEGIN{OFS="\t"} {c[$2]++} END{for (t in c) print t, c[t]}' "${READ2TAXID_TSV}" \
  | LC_ALL=C sort -t $'\t' -k1,1n \
  > "${COUNTS_NUM_TSV}.tmp"

echo -e "taxid\tread_count" > "${COUNTS_NUM_TSV}"
cat "${COUNTS_NUM_TSV}.tmp" >> "${COUNTS_NUM_TSV}"
rm -f "${COUNTS_NUM_TSV}.tmp"

N_NUM_LINES=$(wc -l < "${COUNTS_NUM_TSV}" || echo 0)
echo "  [INFO] Taxa (incl. header): ${N_NUM_LINES}"

if [[ "${N_NUM_LINES}" -le 1 ]]; then
    echo "WARNING: No taxa after filtering. Writing empty taxonomy table."
    echo -e "taxid\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
    exit 0
fi

########################
# 3) TAXID -> LINEAGE + CANONICAL RANKS (TaxonKit)
########################
echo "[$(date)] Building unique taxid list -> ${TAXIDS_TXT}"
cut -f1 "${COUNTS_NUM_TSV}" | tail -n +2 | LC_ALL=C sort -u > "${TAXIDS_TXT}"

echo "[$(date)] Building taxid -> lineage (names) -> ${TAXLINEAGE}"
taxonkit lineage -n "${TAXIDS_TXT}" > "${TAXLINEAGE}.tmp"
echo -e "taxid\tlineage" > "${TAXLINEAGE}"
cat "${TAXLINEAGE}.tmp" >> "${TAXLINEAGE}"
rm -f "${TAXLINEAGE}.tmp"

# Key findings:
# - Never parse lineage by position.
# - Use {d} (domain/superkingdom) instead of {k} (kingdom) because {k} may be blank.
# - Ensure taxid stays the join key by prepending it explicitly.
echo "[$(date)] Building taxid -> canonical ranks (Domain..Species) -> ${TAXRANKS}"
taxonkit reformat -I 1 -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" "${TAXIDS_TXT}" > "${TAXRANKS}.ranks.tmp"

echo -e "taxid\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXRANKS}"
paste "${TAXIDS_TXT}" "${TAXRANKS}.ranks.tmp" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
  >> "${TAXRANKS}"
rm -f "${TAXRANKS}.ranks.tmp"

# Merge lineage + ranks into one map keyed by taxid
LC_ALL=C tail -n +2 "${TAXLINEAGE}" | sort -t $'\t' -k1,1 > "${TAXLINEAGE}.sorted"
LC_ALL=C tail -n +2 "${TAXRANKS}"   | sort -t $'\t' -k1,1 > "${TAXRANKS}.sorted"

LC_ALL=C join -t $'\t' -1 1 -2 1 \
  "${TAXLINEAGE}.sorted" \
  "${TAXRANKS}.sorted" \
  > "${TAXMAP}.tmp"
# taxid  lineage  Domain  phylum  class  order  family  genus  species

echo -e "taxid\tlineage\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXMAP}"
cat "${TAXMAP}.tmp" >> "${TAXMAP}"
rm -f "${TAXMAP}.tmp" "${TAXLINEAGE}.sorted" "${TAXRANKS}.sorted"

########################
# 4) JOIN COUNTS + TAXONOMY -> FINAL TABLE
########################
echo "[$(date)] Building final taxonomy table -> ${FINAL_TABLE}"

LC_ALL=C tail -n +2 "${COUNTS_NUM_TSV}" | sort -t $'\t' -k1,1 > "${COUNTS_NUM_TSV}.sorted"
LC_ALL=C tail -n +2 "${TAXMAP}"         | sort -t $'\t' -k1,1 > "${TAXMAP}.sorted"

LC_ALL=C join -t $'\t' -1 1 -2 1 \
  "${COUNTS_NUM_TSV}.sorted" \
  "${TAXMAP}.sorted" \
  > "${FINAL_TABLE}.joined"
# taxid  read_count  lineage  Domain  phylum  class  order  family  genus  species

echo -e "taxid\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"

awk -F'\t' 'BEGIN{OFS="\t"}
{
    taxid=$1
    if (taxid=="1") next

    count=$2
    lineage=$3
    D=$4; P=$5; C=$6; O=$7; F=$8; G=$9; S=$10

    rank=""
    if (S!="")      rank="species"
    else if (G!="") rank="genus"
    else if (F!="") rank="family"
    else if (O!="") rank="order"
    else if (C!="") rank="class"
    else if (P!="") rank="phylum"
    else if (D!="") rank="domain"

    print taxid, count, lineage, rank, D, P, C, O, F, G, S
}' "${FINAL_TABLE}.joined" >> "${FINAL_TABLE}"

rm -f "${COUNTS_NUM_TSV}.sorted" "${TAXMAP}.sorted" "${FINAL_TABLE}.joined"

echo "[$(date)] Finished taxonomy table for sample: ${SAMPLE}"
echo "Filtered Kraken file : ${KRAKEN_FILTERED}"
echo "Final table          : ${FINAL_TABLE}"

