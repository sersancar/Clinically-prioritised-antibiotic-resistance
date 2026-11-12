#!/bin/bash
#SBATCH --job-name=CARD_ARG_Tax_fromBAM
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=1-60
#SBATCH --output=/home/ssanchez/AREST_scripts/logs/%x_%A_%a.out
#SBATCH --error=/home/ssanchez/AREST_scripts/logs/%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sergio.sanchezcarrillo@universityofgalway.ie

set -euo pipefail

########################
# USER CONFIGURATION
########################

SAMPLE_LIST_FILE="/home/ssanchez/AREST_scripts/sample_names.txt"

# Existing CARD mapping results
CARD_DIR="/data3/ssanchez/card_genes"
BAM_DIR="${CARD_DIR}/alignments"
BAM_SUFFIX=".filt.sorted.bam"

# CARD reference (only needed for fallback coverage calculation)
REFERENCEFILE="${CARD_DIR}/card_genes.fna"

# Pre-filtered Kraken2 outputs (already score-filtered at 0.2)
KRAKEN_FILTERED_DIR="/data3/ssanchez/ARG_Taxon_link"
KRAKEN_FILTERED_SUFFIX=".kraken.score0.2.out"
KRAKEN_SCORE=0.2

# Output directory for linkage tables
OUTDIR="${CARD_DIR}/CARD_ARG_Taxon_link"
mkdir -p "${OUTDIR}"

# Evidence thresholds for calling ARG “present”
MIN_MAPQ=30
MIN_BREADTH=0.80
MIN_MEAN_DEPTH=5

# Taxonkit DB:
# IMPORTANT: point to the directory that contains nodes.dmp and names.dmp.
# In your case, that is the Kraken DB root.
export TAXONKIT_DB="/data3/ssanchez/kraken2DB"

# Threads
SAMTOOLS_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Wait until main disk is mounted (optional safety)
while [ ! -d "/data3/ssanchez" ]; do
    echo "Directory /data3/ssanchez does not exist. Waiting..."
    sleep 5
done

########################
# Conda env (samtools, taxonkit)
########################
source /data3/ssanchez/miniforge3/etc/profile.d/conda.sh
# conda activate <your_env>

command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not in PATH" >&2; exit 1; }
command -v taxonkit  >/dev/null 2>&1 || { echo "ERROR: taxonkit not in PATH" >&2; exit 1; }

echo "Job started on $(date)"
echo "Node(s): ${SLURM_NODELIST:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "BAM dir: ${BAM_DIR}"
echo "Kraken dir (filtered): ${KRAKEN_FILTERED_DIR} (score=${KRAKEN_SCORE})"
echo "Thresholds: MAPQ>=${MIN_MAPQ}, breadth>=${MIN_BREADTH}, mean_depth>=${MIN_MEAN_DEPTH}"
echo "Output dir: ${OUTDIR}"
echo "TAXONKIT_DB: ${TAXONKIT_DB}"

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

########################
# INPUTS
########################
BAM="${BAM_DIR}/${SAMPLE}${BAM_SUFFIX}"
KRAKEN_FILTERED="${KRAKEN_FILTERED_DIR}/${SAMPLE}${KRAKEN_FILTERED_SUFFIX}"

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}" >&2
    exit 1
fi
if [[ ! -f "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found: ${BAM}.bai" >&2
    exit 1
fi
if [[ ! -f "${KRAKEN_FILTERED}" ]]; then
    echo "ERROR: Filtered Kraken file not found: ${KRAKEN_FILTERED}" >&2
    exit 1
fi

# For fallback coverage mode only
if [[ -f "${REFERENCEFILE}" ]]; then
    if [[ ! -f "${REFERENCEFILE}.fai" ]]; then
        samtools faidx "${REFERENCEFILE}"
    fi
else
    echo "WARNING: Reference FASTA not found (${REFERENCEFILE}). Fallback coverage mode will not be available." >&2
fi

########################
# OUTPUTS (per sample)
########################
PREFIX="${OUTDIR}/${SAMPLE}"

ARG_METRICS="${PREFIX}.ARG_metrics.tsv"
ARG_PRESENT="${PREFIX}.ARG_present.tsv"
PRESENT_GENES="${PREFIX}.present_genes.txt"

READ2ARG_TSV="${PREFIX}.read_to_arg.tsv"
READ2TAXID_TSV="${PREFIX}.read_to_taxid.tsv"
READ_ARG_TAXID_TSV="${PREFIX}.read_arg_taxid.tsv"

COUNTS_NUM_TSV="${PREFIX}.CARD_ARG_taxid_counts.numeric.tsv"

# Taxonomy maps (audit + correct ranks)
TAXIDS_TXT="${PREFIX}.taxids.txt"
TAXLINEAGE="${PREFIX}.taxid_lineage.tsv"   # taxid + full lineage (audit)
TAXRANKS="${PREFIX}.taxid_ranks.tsv"       # taxid + Domain..Species
TAXMAP="${PREFIX}.taxid_taxonomy.tsv"      # taxid + lineage + Domain..Species

FINAL_TABLE="${PREFIX}.CARD_ARG_taxonomy.tsv"

########################
# 1) PER-GENE METRICS + PRESENT LIST (from existing BAM)
########################
echo -e "gene\tlength\tcov_bases\tbreadth\tmean_depth\tnum_reads" > "${ARG_METRICS}"
echo -e "gene\tlength\tcov_bases\tbreadth\tmean_depth\tnum_reads" > "${ARG_PRESENT}"

if samtools coverage -h >/dev/null 2>&1; then
  echo "[$(date)] Computing per-gene coverage with samtools coverage..."
  samtools coverage "${BAM}" \
    | awk -v OFS="\t" -v mb="${MIN_BREADTH}" -v md="${MIN_MEAN_DEPTH}" '
      $0 ~ /^#/ {next}
      NF < 7 {next}
      {
        gene=$1; start=$2; end=$3; numreads=$4; covbases=$5; meandepth=$7;
        len=end-start+1;
        breadth=(len>0 ? covbases/len : 0);
        print gene,len,covbases,breadth,meandepth,numreads >> "'"${ARG_METRICS}"'"
        if (breadth>=mb && meandepth>=md) {
          print gene,len,covbases,breadth,meandepth,numreads >> "'"${ARG_PRESENT}"'"
        }
      }'
else
  if [[ ! -f "${REFERENCEFILE}.fai" ]]; then
    echo "ERROR: samtools coverage not available and ${REFERENCEFILE}.fai not present for fallback." >&2
    exit 1
  fi
  echo "[$(date)] samtools coverage not available; falling back to samtools depth (-aa) + .fai"
  samtools depth -aa -d 0 "${BAM}" \
    | awk -v OFS="\t" -v fai="${REFERENCEFILE}.fai" -v mb="${MIN_BREADTH}" -v md="${MIN_MEAN_DEPTH}" '
        BEGIN { while ((getline < fai) > 0) { len[$1]=$2 } close(fai) }
        { g=$1; d=$3; sumd[g]+=d; if (d>0) covpos[g]++ }
        END {
          for (g in len) {
            L=len[g]
            cb=(covpos[g]+0)
            breadth=(L>0 ? cb/L : 0)
            meandepth=(L>0 ? sumd[g]/L : 0)
            print g,L,cb,breadth,meandepth,"NA" >> "'"${ARG_METRICS}"'"
            if (breadth>=mb && meandepth>=md) {
              print g,L,cb,breadth,meandepth,"NA" >> "'"${ARG_PRESENT}"'"
            }
          }
        }'
fi

tail -n +2 "${ARG_PRESENT}" | cut -f1 | LC_ALL=C sort -u > "${PRESENT_GENES}"
N_PRESENT=$(wc -l < "${PRESENT_GENES}" || true)
echo "[$(date)] Present genes passing thresholds: ${N_PRESENT}"

if [[ "${N_PRESENT}" -eq 0 ]]; then
  echo "WARNING: No ARGs pass presence thresholds for ${SAMPLE}. Writing empty final table."
  echo -e "taxid\tARG_ID\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
  exit 0
fi

########################
# 2) READ -> ARG (from existing BAM)
########################
echo "[$(date)] Extracting read-to-ARG mapping (MAPQ>=${MIN_MAPQ}) -> ${READ2ARG_TSV}"

samtools view -@ "${SAMTOOLS_THREADS}" -F 2308 -q "${MIN_MAPQ}" "${BAM}" \
  | awk 'BEGIN{OFS="\t"} {print $1,$3}' \
  | LC_ALL=C sort -u \
  > "${READ2ARG_TSV}"
# cols: readID  ARG_ID

########################
# 3) READ -> TAXID (from pre-filtered Kraken file)
########################
echo "[$(date)] Extracting read-to-taxid from ${KRAKEN_FILTERED} -> ${READ2TAXID_TSV}"

awk 'BEGIN{FS="\t"; OFS="\t"}
     $1=="C" {
        rid=$2; f=$3;
        taxid=""
        if (f ~ /^[0-9]+$/) taxid=f
        else if (match(f, /\(taxid ([0-9]+)\)/, m)) taxid=m[1]
        else next
        if (taxid=="0" || taxid=="1") next
        print rid, taxid
     }' "${KRAKEN_FILTERED}" \
  | LC_ALL=C sort -u \
  > "${READ2TAXID_TSV}"
# cols: readID  taxid

########################
# 4) JOIN READ->ARG + READ->TAXID
########################
echo "[$(date)] Joining read_to_arg and read_to_taxid -> ${READ_ARG_TAXID_TSV}"

LC_ALL=C join -t $'\t' -1 1 -2 1 \
  <(LC_ALL=C sort -t $'\t' -k1,1 "${READ2ARG_TSV}") \
  <(LC_ALL=C sort -t $'\t' -k1,1 "${READ2TAXID_TSV}") \
  > "${READ_ARG_TAXID_TSV}"
# cols: readID  ARG_ID  taxid

if [[ ! -s "${READ_ARG_TAXID_TSV}" ]]; then
  echo "WARNING: No read-level overlaps between BAM and Kraken for ${SAMPLE}."
  echo -e "taxid\tARG_ID\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
  exit 0
fi

########################
# 5) AGGREGATE ARG x TAXID COUNTS (ONLY PRESENT GENES)
########################
echo "[$(date)] Aggregating ARG x taxid read counts (only present genes) -> ${COUNTS_NUM_TSV}"

awk -v OFS="\t" '
NR==FNR { present[$1]=1; next }
{
  arg=$2; tax=$3
  if (!present[arg]) next
  key=arg "\t" tax
  c[key]++
}
END { for (k in c) print k, c[k] }
' "${PRESENT_GENES}" "${READ_ARG_TAXID_TSV}" \
  | LC_ALL=C sort -t $'\t' -k1,1 -k2,2 \
  > "${COUNTS_NUM_TSV}.tmp"

echo -e "ARG_ID\ttaxid\tread_count" > "${COUNTS_NUM_TSV}"
cat "${COUNTS_NUM_TSV}.tmp" >> "${COUNTS_NUM_TSV}"
rm -f "${COUNTS_NUM_TSV}.tmp"

N_NUM_LINES=$(wc -l < "${COUNTS_NUM_TSV}" || true)
echo "  [INFO] ARG–taxid rows (incl. header): ${N_NUM_LINES}"

if [[ "${N_NUM_LINES}" -le 1 ]]; then
  echo "WARNING: No ARG–taxid pairs after filtering to present genes for ${SAMPLE}."
  echo -e "taxid\tARG_ID\tread_count\tlineage\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TABLE}"
  exit 0
fi

########################
# 6) TAXID -> TAXONOMY (CORRECT WAY)
#    - Full lineage for auditing: taxonkit lineage -n
#    - Canonical ranks: taxonkit reformat -I 1
#      IMPORTANT: prepend taxid explicitly; do NOT parse by position.
#      Use {d} (domain/superkingdom) rather than {k} (kingdom) because {k} is often blank.
########################
echo "[$(date)] Building taxid list -> ${TAXIDS_TXT}"
cut -f2 "${COUNTS_NUM_TSV}" | tail -n +2 | LC_ALL=C sort -u > "${TAXIDS_TXT}"

echo "[$(date)] Building taxid -> lineage (names) -> ${TAXLINEAGE}"
taxonkit lineage -n "${TAXIDS_TXT}" > "${TAXLINEAGE}.tmp"
echo -e "taxid\tlineage" > "${TAXLINEAGE}"
cat "${TAXLINEAGE}.tmp" >> "${TAXLINEAGE}"
rm -f "${TAXLINEAGE}.tmp"

echo "[$(date)] Building taxid -> canonical ranks (Domain..Species) -> ${TAXRANKS}"
taxonkit reformat -I 1 -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" "${TAXIDS_TXT}" > "${TAXRANKS}.ranks.tmp"

echo -e "taxid\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXRANKS}"
paste "${TAXIDS_TXT}" "${TAXRANKS}.ranks.tmp" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
  >> "${TAXRANKS}"
rm -f "${TAXRANKS}.ranks.tmp"

# Combine lineage + ranks into one map
LC_ALL=C tail -n +2 "${TAXLINEAGE}" | sort -t $'\t' -k1,1 > "${TAXLINEAGE}.sorted"
LC_ALL=C tail -n +2 "${TAXRANKS}"   | sort -t $'\t' -k1,1 > "${TAXRANKS}.sorted"

LC_ALL=C join -t $'\t' -1 1 -2 1 \
  "${TAXLINEAGE}.sorted" \
  "${TAXRANKS}.sorted" \
  > "${TAXMAP}.tmp"
# join output: taxid  lineage  Domain  phylum  class  order  family  genus  species

echo -e "taxid\tlineage\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXMAP}"
cat "${TAXMAP}.tmp" >> "${TAXMAP}"
rm -f "${TAXMAP}.tmp" "${TAXLINEAGE}.sorted" "${TAXRANKS}.sorted"

########################
# 7) JOIN COUNTS + TAXMAP -> FINAL TABLE
########################
echo "[$(date)] Joining ARG-taxid counts with taxonomy -> ${FINAL_TABLE}"

LC_ALL=C tail -n +2 "${COUNTS_NUM_TSV}" | sort -t $'\t' -k2,2 > "${COUNTS_NUM_TSV}.sorted"
LC_ALL=C tail -n +2 "${TAXMAP}"         | sort -t $'\t' -k1,1 > "${TAXMAP}.sorted"

LC_ALL=C join -t $'\t' -1 2 -2 1 \
  "${COUNTS_NUM_TSV}.sorted" \
  "${TAXMAP}.sorted" \
  > "${FINAL_TABLE}.joined"
# join output:
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
echo "  Final ARG-taxonomy table: ${FINAL_TABLE}"

