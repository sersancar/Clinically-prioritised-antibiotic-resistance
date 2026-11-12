#!/bin/bash
#SBATCH --job-name=AMR_matrix_MAPQ30
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --array=1-42
#SBATCH --output=/home/ncahill/sergio/logs/%x_%A_%a.out
#SBATCH --error=/home/ncahill/sergio/logs/%x_%A_%a.err

set -euo pipefail

########################
# USER CONFIGURATION
########################

SAMPLE_LIST_FILE="/home/ncahill/sergio/samplesBatch1.txt"

AMR_BAM_DIR="/data4/ncahill/AMR++Out/Alignment/BAM_files/Deduped"
BAM_SUFFIX="_alignment_dedup.bam"

MIN_MAPQ=30
RESISTOME_THRESHOLD=80
SAMTOOLS_THREADS=8

# Exclude:
# 4    unmapped
# 256  secondary
# 1024 duplicate
# 2048 supplementary
EXCLUDE_FLAGS=3332

OUTDIR="/data4/ncahill/AMR++Out/filteredResults"
mkdir -p "${OUTDIR}"

MATRIX_CSV="${OUTDIR}/dedup_AMR_analytic_matrix.csv"
GENE_LIST="${OUTDIR}/gene_list.mapq${MIN_MAPQ}.txt"

source /data4/ncahill/miniconda3/etc/profile.d/conda.sh
# conda activate <your_env_with_samtools>

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not in PATH" >&2
    exit 1
}

if [[ ! -f "${SAMPLE_LIST_FILE}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

NSAMPLES=$(grep -cvE '^\s*($|#)' "${SAMPLE_LIST_FILE}")

if [[ "${NSAMPLES}" -le 0 ]]; then
    echo "ERROR: No valid sample names found in ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

echo "Job started on $(date)"
echo "Running on node(s): ${SLURM_NODELIST:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Using BAM dir           : ${AMR_BAM_DIR}"
echo "MAPQ cutoff             : ${MIN_MAPQ}"
echo "Resistome threshold     : ${RESISTOME_THRESHOLD}%"
echo "Exclude flags           : ${EXCLUDE_FLAGS}"
echo "Output dir              : ${OUTDIR}"
echo "Sample list file        : ${SAMPLE_LIST_FILE}"
echo "NSAMPLES inferred       : ${NSAMPLES}"

########################
# GET SAMPLE FOR THIS ARRAY TASK
########################

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit with:" >&2
    echo "  sbatch --array=1-${NSAMPLES} $(basename "$0")" >&2
    exit 1
fi

LINE_NO="${SLURM_ARRAY_TASK_ID}"

SAMPLE="$(awk 'NF && $0 !~ /^#/ {++n; if (n=='"${LINE_NO}"') {print; exit}}' "${SAMPLE_LIST_FILE}")"

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample name for array index ${LINE_NO}" >&2
    exit 1
fi

echo "Processing sample: ${SAMPLE}"

########################
# PER-SAMPLE COVERAGE + GENE COUNTS
########################

BAM="${AMR_BAM_DIR}/${SAMPLE}${BAM_SUFFIX}"

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found for sample ${SAMPLE}: ${BAM}" >&2
    exit 1
fi

COV_FILE="${OUTDIR}/${SAMPLE}.mapq${MIN_MAPQ}.coverage.tsv"
SAMPLE_COUNTS="${OUTDIR}/${SAMPLE}.mapq${MIN_MAPQ}.gene_counts.tsv"
TMP_COUNTS="${SAMPLE_COUNTS}.tmp"

echo "  BAM: ${BAM}"
echo "  Coverage file: ${COV_FILE}"
echo "  Output counts: ${SAMPLE_COUNTS}"

########## 1) Compute coverage fraction per gene ##########

samtools depth \
    -@ "${SAMTOOLS_THREADS}" \
    -aa \
    -Q "${MIN_MAPQ}" \
    -G "${EXCLUDE_FLAGS}" \
    "${BAM}" \
    | awk 'BEGIN{OFS="\t"}
           {
               gene = $1
               depth = $3
               len[gene]++
               if (depth > 0) cov[gene]++
           }
           END {
               for (g in len) {
                   frac = (cov[g] / len[g]) * 100
                   printf "%s\t%.2f\n", g, frac
               }
           }' \
    > "${COV_FILE}"

########## 2) Count reads per gene and apply coverage filter ##########

echo -e "gene_accession\tcount" > "${TMP_COUNTS}"

samtools view \
    -@ "${SAMTOOLS_THREADS}" \
    -F "${EXCLUDE_FLAGS}" \
    -q "${MIN_MAPQ}" \
    "${BAM}" \
    | awk '{print $3}' \
    | LC_ALL=C sort \
    | uniq -c \
    | awk -v thresh="${RESISTOME_THRESHOLD}" '
        NR==FNR {
            cov[$1] = $2
            next
        }
        {
            count = $1
            gene  = $2
            if (cov[gene] + 0 >= thresh)
                printf "%s\t%d\n", gene, count
        }' "${COV_FILE}" - \
    >> "${TMP_COUNTS}"

mv "${TMP_COUNTS}" "${SAMPLE_COUNTS}"

echo "Finished counting for sample ${SAMPLE} at $(date)"

########################
# MERGE LOGIC: ONLY LAST TASK DOES IT
########################

if [[ "${SLURM_ARRAY_TASK_ID}" -ne "${NSAMPLES}" ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID}: counts done. Merge will be done by task ${NSAMPLES}."
    exit 0
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: will perform final merge after all counts are ready."

declare -a EXPECTED_COUNTS

while read -r s; do
    [[ -z "${s}" ]] && continue
    [[ "${s}" =~ ^# ]] && continue
    EXPECTED_COUNTS+=( "${OUTDIR}/${s}.mapq${MIN_MAPQ}.gene_counts.tsv" )
done < "${SAMPLE_LIST_FILE}"

echo "Waiting for all ${#EXPECTED_COUNTS[@]} count files to appear..."

while true; do
    missing=0

    for f in "${EXPECTED_COUNTS[@]}"; do
        if [[ ! -s "${f}" ]]; then
            missing=1
            break
        fi
    done

    if [[ "${missing}" -eq 0 ]]; then
        echo "All count files found at $(date). Proceeding to merge."
        break
    fi

    echo "Not all count files present yet. Sleeping 60 seconds..."
    sleep 60
done

########################
# COLLECT COUNT FILES IN SAMPLE ORDER
########################

SAMPLES=()
COUNT_FILES=()

while read -r s; do
    [[ -z "${s}" ]] && continue
    [[ "${s}" =~ ^# ]] && continue

    counts="${OUTDIR}/${s}.mapq${MIN_MAPQ}.gene_counts.tsv"

    if [[ -f "${counts}" ]]; then
        SAMPLES+=( "${s}" )
        COUNT_FILES+=( "${counts}" )
    else
        echo "WARNING: counts file not found for sample ${s}: ${counts}" >&2
    fi
done < "${SAMPLE_LIST_FILE}"

if [[ ${#COUNT_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No per-sample count files found in ${OUTDIR}. Aborting merge." >&2
    exit 1
fi

echo "Samples to include in matrix: ${#SAMPLES[@]}"
printf '  - %s\n' "${SAMPLES[@]}"

########################
# BUILD UNION OF GENES
########################

TMP_ALL="${OUTDIR}/all_genes.mapq${MIN_MAPQ}.tmp"
: > "${TMP_ALL}"

echo "[$(date)] Building union of AMR gene_accessions -> ${GENE_LIST}"

for f in "${COUNT_FILES[@]}"; do
    if [[ -s "${f}" ]]; then
        tail -n +2 "${f}" >> "${TMP_ALL}"
    fi
done

if [[ ! -s "${TMP_ALL}" ]]; then
    echo "WARNING: No genes found. Writing empty matrix."
    echo "gene_accession" > "${MATRIX_CSV}"
    rm -f "${TMP_ALL}"
    exit 0
fi

cut -f1 "${TMP_ALL}" | LC_ALL=C sort -u > "${GENE_LIST}"
rm -f "${TMP_ALL}"

N_GENES=$(wc -l < "${GENE_LIST}" || echo 0)

echo "  [INFO] Total unique AMR gene_accessions: ${N_GENES}"

if [[ "${N_GENES}" -le 0 ]]; then
    echo "WARNING: No genes found after union. Writing empty matrix."
    echo "gene_accession" > "${MATRIX_CSV}"
    exit 0
fi

########################
# MERGE INTO WIDE MATRIX
########################

echo "[$(date)] Merging per-sample counts into matrix -> ${MATRIX_CSV}"

SAMPLE_LIST_STR=""

for s in "${SAMPLES[@]}"; do
    SAMPLE_LIST_STR+="${s} "
done

awk -v samples="${SAMPLE_LIST_STR}" '
BEGIN {
    FS = OFS = "\t"
    split(samples, sname, " ")
    file_idx = 0
}

NR == FNR {
    genes[++ng] = $1
    next
}

FNR == 1 && NR != 1 {
    file_idx++
    next
}

{
    gene = $1
    count = $2 + 0
    acc[gene, file_idx] = count
}

END {
    sep = ","
    nfiles = file_idx

    printf "gene_accession"

    for (i = 1; i <= nfiles; i++) {
        if (sname[i] == "") sname[i] = "sample" i
        printf "%s%s", sep, sname[i]
    }

    printf "\n"

    for (g = 1; g <= ng; g++) {
        gene = genes[g]
        printf "%s", gene

        for (i = 1; i <= nfiles; i++) {
            key = gene SUBSEP i
            val = (key in acc ? acc[key] : 0)
            printf "%s%s", sep, val
        }

        printf "\n"
    }
}
' "${GENE_LIST}" "${COUNT_FILES[@]}" > "${MATRIX_CSV}"

echo "[$(date)] Finished merge."
echo "Final matrix: ${MATRIX_CSV}"
echo "Dimensions: ${N_GENES} genes x ${#SAMPLES[@]} samples plus header"
