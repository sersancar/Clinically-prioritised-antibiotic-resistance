#!/bin/bash
#SBATCH --job-name=ARG_Tax_linkage
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=1-60
#SBATCH --output=/home/ssanchez/AREST_scripts/logs/%x_%A_%a.out
#SBATCH --error=/home/ssanchez/AREST_scripts/logs/%x_%A_%a.err

set -euo pipefail

########################
# USER CONFIGURATION
########################

# File with sample names (one per line, no header), 60 lines
SAMPLE_LIST_FILE="/home/ssanchez/AREST_scripts/sample_names.txt"

# AMR++ deduplicated BAM directory
AMR_BAM_DIR="/data3/ssanchez/AMR++Results/AMR++Out/Alignment/BAM_files/Deduped"
BAM_SUFFIX="_alignment_dedup.bam"

# Kraken2 results directory (UNFILTERED .out; this script filters by score)
KRAKEN_RES_DIR="/data3/ssanchez/kraken2"
KRAKEN_SUFFIX=".out"

# Output directory
OUTDIR="/data3/ssanchez/ARG_Taxon_link"
mkdir -p "${OUTDIR}"

# Control strength of ARG alignment (mapping quality threshold)
MIN_MAPQ=30

# AMR++ resistome threshold (minimum % of gene length covered by depth>0)
RESISTOME_THRESHOLD=80

# Control strength of Kraken2 taxonomic assignments (score on k-mer field)
KRAKEN_SCORE=0.2

# TaxonKit DB
# IMPORTANT: set to the directory that contains nodes.dmp and names.dmp.
# If you are using the Kraken DB layout like /data3/ssanchez/kraken2DB with nodes.dmp/names.dmp in the root,
# set TAXONKIT_DB to that DB root instead of /data3/ssanchez/NCBIDBs.
NCBI_TAXDIR="/data3/ssanchez/kraken2DB"
export TAXONKIT_DB="${NCBI_TAXDIR}"

# Threads for samtools
SAMTOOLS_THREADS="${SLURM_CPUS_PER_TASK:-8}"

########################
# WAIT FOR DISK
########################
while [ ! -d "/data3/ssanchez" ]; do
    echo "Directory /data3/ssanchez does not exist. Waiting..."
    sleep 5
done

########################
# ENV
########################
source /data3/ssanchez/miniforge3/etc/profile.d/conda.sh
# conda activate <your_env_with_samtools_taxonkit_python>

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not in PATH" >&2; exit 1; }
command -v taxonkit >/dev/null 2>&1 || { echo "ERROR: taxonkit not in PATH" >&2; exit 1; }
command -v python   >/dev/null 2>&1 || { echo "ERROR: python not in PATH" >&2; exit 1; }

echo "Job started on $(date)"
echo "Running on node(s): $SLURM_NODELIST"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Using BAM dir          : ${AMR_BAM_DIR}"
echo "Using Kraken dir       : ${KRAKEN_RES_DIR}"
echo "MAPQ cutoff            : ${MIN_MAPQ}"
echo "Resistome threshold    : ${RESISTOME_THRESHOLD}%"
echo "Kraken score cutoff    : ${KRAKEN_SCORE}"
echo "Output dir             : ${OUTDIR}"
echo "Samtools threads       : ${SAMTOOLS_THREADS}"
echo "TAXONKIT_DB            : ${TAXONKIT_DB}"

########################
# GET SAMPLE NAME (ARRAY MODE)
########################
if [[ ! -f "${SAMPLE_LIST_FILE}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST_FILE}" | tr -d '\r' | xargs)"
if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample name on line ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

echo "Processing sample: ${SAMPLE}"

BAM="${AMR_BAM_DIR}/${SAMPLE}${BAM_SUFFIX}"
KRAKEN_ORIG="${KRAKEN_RES_DIR}/${SAMPLE}${KRAKEN_SUFFIX}"

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}" >&2
    exit 1
fi
if [[ ! -f "${KRAKEN_ORIG}" ]]; then
    echo "ERROR: Kraken file not found: ${KRAKEN_ORIG}" >&2
    exit 1
fi

PREFIX="${OUTDIR}/${SAMPLE}"

# Outputs
READ2ARG_TSV="${PREFIX}.read_to_arg.tsv"
COV_FILE="${PREFIX}.mapq${MIN_MAPQ}.coverage.tsv"

KRAKEN_FILTERED="${PREFIX}.kraken.score${KRAKEN_SCORE}.out"
READ2TAXID_TSV="${PREFIX}.read_to_taxid.tsv"

READ_ARG_TAXID_TSV="${PREFIX}.read_arg_taxid.tsv"

COUNTS_NUM_TSV="${PREFIX}.ARG_taxid_counts.numeric.tsv"
TAXIDS_TXT="${PREFIX}.taxids.txt"
TAXLINEAGE="${PREFIX}.taxid_lineage.tsv"
TAXRANKS="${PREFIX}.taxid_ranks.tsv"
TAXMAP="${PREFIX}.taxid_taxonomy.tsv"

FINAL_TABLE="${OUTDIR}/${SAMPLE}.ARG_taxonomy.tsv"

########################
# 1) READ -> ARG (FROM AMR BAM, WITH MAPQ FILTER)
########################
echo "[$(date)] Extracting read-to-ARG mapping with MAPQ >= ${MIN_MAPQ} -> ${READ2ARG_TSV}"

samtools view -@ "${SAMTOOLS_THREADS}" -F 2308 -q "${MIN_MAPQ}" "${BAM}" \
  | awk 'BEGIN{OFS="\t"} {print $1,$3}' \
  | sort -u \
  > "${READ2ARG_TSV}"
# cols: readID  ARG_ID

########################
# 1b) PER-GENE COVERAGE FRACTION (RESISTOME THRESHOLD)
########################
echo "[$(date)] Computing per-gene coverage fraction (MAPQ >= ${MIN_MAPQ}) -> ${COV_FILE}"

samtools depth -aa -@ "${SAMTOOLS_THREADS}" -G 2048 -Q "${MIN_MAPQ}" "${BAM}" \
  | awk 'BEGIN{OFS="\t"}
         {
             gene = $1; depth = $3;
             len[gene]++;               # total positions for this gene
             if (depth > 0) cov[gene]++ # positions with coverage
         }
         END {
             for (g in len) {
                 frac = (cov[g] / len[g]) * 100;
                 printf "%s\t%.2f\n", g, frac;
             }
         }' \
  | sort -k1,1 \
  > "${COV_FILE}"
# cols: ARG_ID   coverage_percent

########################
# 2) FILTER KRAKEN BY SCORE + EXTRACT READ->TAXID (NUMERIC)
########################
echo "[$(date)] Filtering Kraken2 by score >= ${KRAKEN_SCORE} -> ${KRAKEN_FILTERED}"
echo "[$(date)] Building read->taxid (numeric) from filtered Kraken2 -> ${READ2TAXID_TSV}"

python - "$KRAKEN_SCORE" "${KRAKEN_ORIG}" "${KRAKEN_FILTERED}" "${READ2TAXID_TSV}" << 'EOF'
import sys, re

min_score = float(sys.argv[1])
infile = sys.argv[2]
filtered_out = sys.argv[3]
read2taxid_out = sys.argv[4]

taxid_pattern = re.compile(r"\(taxid\s+(\d+)\)")

def extract_label_taxid(field: str) -> str:
    s = field.strip()
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
            continue

        read_id = parts[1]
        label_field = parts[2]
        label_taxid = extract_label_taxid(label_field)
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

########################
# 3) JOIN READ->ARG + READ->TAXID
########################
echo "[$(date)] Joining read_to_arg and read_to_taxid -> ${READ_ARG_TAXID_TSV}"

LC_ALL=C join -t $'\t' -1 1 -2 1 \
  <(LC_ALL=C sort -t $'\t' -k1,1 "${READ2ARG_TSV}") \
  <(LC_ALL=C sort -t $'\t' -k1,1 "${READ2TAXID_TSV}") \
  > "${READ_ARG_TAXID_TSV}"
# cols: readID  ARG_ID  taxid

if [[ ! -s "${READ_ARG_TAXID_TSV}" ]]; then
    echo "WARNING: No read-level overlaps between AMR BAM and Kraken after filtering for ${SAMPLE}."
    echo -e "taxid\tARG_ID\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
    echo "[$(date)] Finished sample (no overlaps): ${SAMPLE}"
    exit 0
fi

########################
# 4) AGGREGATE ARG x TAXID COUNTS, APPLY COVERAGE THRESHOLD
########################
echo "[$(date)] Aggregating ARG x taxid counts with coverage >= ${RESISTOME_THRESHOLD}% -> ${COUNTS_NUM_TSV}"

awk -v thresh="${RESISTOME_THRESHOLD}" -F'\t' '
NR==FNR { cov[$1]=$2; next }
{
  arg=$2; tax=$3
  if ((cov[arg] + 0) < thresh) next
  key=arg "\t" tax
  c[key]++
}
END {
  for (k in c) print k "\t" c[k]
}
' "${COV_FILE}" "${READ_ARG_TAXID_TSV}" \
  | LC_ALL=C sort -t $'\t' -k1,1 -k2,2 \
  > "${COUNTS_NUM_TSV}.tmp"

echo -e "ARG_ID\ttaxid\tread_count" > "${COUNTS_NUM_TSV}"
cat "${COUNTS_NUM_TSV}.tmp" >> "${COUNTS_NUM_TSV}"
rm -f "${COUNTS_NUM_TSV}.tmp"

N_NUM_LINES=$(wc -l < "${COUNTS_NUM_TSV}" || true)
echo "  [INFO] Numeric ARG–taxid rows (incl. header): ${N_NUM_LINES}"

if [[ "${N_NUM_LINES}" -le 1 ]]; then
    echo "WARNING: No ARG–taxid pairs after coverage filter. Writing empty final table."
    echo -e "taxid\tARG_ID\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
    echo "[$(date)] Finished sample (empty counts): ${SAMPLE}"
    exit 0
fi

########################
# 5) TAXID -> TAXONOMY (Domain..Species; explicit taxid column)
########################
echo "[$(date)] Building taxid list -> ${TAXIDS_TXT}"
cut -f2 "${COUNTS_NUM_TSV}" | tail -n +2 | sort -u > "${TAXIDS_TXT}"

echo "[$(date)] Building taxid -> lineage (names) -> ${TAXLINEAGE}"
taxonkit lineage -n "${TAXIDS_TXT}" > "${TAXLINEAGE}.tmp"
echo -e "taxid\tlineage" > "${TAXLINEAGE}"
cat "${TAXLINEAGE}.tmp" >> "${TAXLINEAGE}"
rm -f "${TAXLINEAGE}.tmp"

# Use {d} (superkingdom/domain) instead of {k} (kingdom) because {k} is often blank in NCBI
echo "[$(date)] Building taxid -> canonical ranks (Domain..Species) -> ${TAXRANKS}"
taxonkit reformat -I 1 -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" "${TAXIDS_TXT}" > "${TAXRANKS}.ranks.tmp"

echo -e "taxid\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXRANKS}"
paste "${TAXIDS_TXT}" "${TAXRANKS}.ranks.tmp" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
  >> "${TAXRANKS}"
rm -f "${TAXRANKS}.ranks.tmp"

# Merge lineage + ranks into a single taxonomy map keyed by taxid
LC_ALL=C tail -n +2 "${TAXLINEAGE}" | sort -t $'\t' -k1,1 > "${TAXLINEAGE}.sorted"
LC_ALL=C tail -n +2 "${TAXRANKS}"   | sort -t $'\t' -k1,1 > "${TAXRANKS}.sorted"

LC_ALL=C join -t $'\t' -1 1 -2 1 \
  "${TAXLINEAGE}.sorted" \
  "${TAXRANKS}.sorted" \
  > "${TAXMAP}.tmp"
# join output:
# taxid  lineage  Domain  phylum  class  order  family  genus  species

echo -e "taxid\tlineage\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXMAP}"
cat "${TAXMAP}.tmp" >> "${TAXMAP}"
rm -f "${TAXMAP}.tmp" "${TAXLINEAGE}.sorted" "${TAXRANKS}.sorted"

########################
# 6) JOIN COUNTS + TAXONOMY -> FINAL TABLE
########################
echo "[$(date)] Joining counts with taxonomy -> ${FINAL_TABLE}"

LC_ALL=C tail -n +2 "${COUNTS_NUM_TSV}" | sort -t $'\t' -k2,2 > "${COUNTS_NUM_TSV}.sorted"
LC_ALL=C tail -n +2 "${TAXMAP}"         | sort -t $'\t' -k1,1 > "${TAXMAP}.sorted"

LC_ALL=C join -t $'\t' -1 2 -2 1 \
  "${COUNTS_NUM_TSV}.sorted" \
  "${TAXMAP}.sorted" \
  > "${FINAL_TABLE}.joined"
# join output columns:
# taxid  ARG_ID  read_count  lineage  Domain  phylum  class  order  family  genus  species

echo -e "taxid\tARG_ID\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"

awk -F'\t' 'BEGIN{OFS="\t"}
{
  taxid=$1
  if (taxid=="1") next

  argid=$2
  count=$3
  lineage=$4
  D=$5; P=$6; C=$7; O=$8; F=$9; G=$10; S=$11

  rank=""
  if (S!="")      rank="species"
  else if (G!="") rank="genus"
  else if (F!="") rank="family"
  else if (O!="") rank="order"
  else if (C!="") rank="class"
  else if (P!="") rank="phylum"
  else if (D!="") rank="domain"

  print taxid,argid,count,lineage,rank,D,P,C,O,F,G,S
}' "${FINAL_TABLE}.joined" >> "${FINAL_TABLE}"

rm -f "${COUNTS_NUM_TSV}.sorted" "${TAXMAP}.sorted" "${FINAL_TABLE}.joined"

echo "[$(date)] Finished sample: ${SAMPLE}"

