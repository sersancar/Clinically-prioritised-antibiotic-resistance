#!/bin/bash
#SBATCH --job-name=AMR_matrix_MAPQ30
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --array=1-60
#SBATCH --output=logs/AMR_matrix_MAPQ30_%A_%a.out
#SBATCH --error=logs/AMR_matrix_MAPQ30_%A_%a.err

set -euo pipefail

########################
# USER CONFIGURATION
########################

# File with sample names (one per line, no header), 60 lines
SAMPLE_LIST_FILE="/home/ssanchez/AREST_scripts/sample_names.txt"

# AMR++ deduplicated BAM directory
AMR_BAM_DIR="/data3/ssanchez/AMR++Results/AMR++Out/Alignment/BAM_files/Deduped"

# BAM filename pattern
BAM_SUFFIX="_alignment_dedup.bam"

# Mapping quality threshold (same as ARG_Tax_linkage.sh)
MIN_MAPQ=30

# AMR++ resistome threshold: minimum % of gene length covered (depth > 0)
RESISTOME_THRESHOLD=80

# Wait until the disk is mounted
while [ ! -d "/data3/ssanchez" ]; do
    echo "Directory does not exist. Waiting..."
    sleep 5
done

# Output directory for matrix and intermediate files
OUTDIR="/data3/ssanchez/AMR++Results/AMR++Out/Results/filteredResults"
mkdir -p "${OUTDIR}"

MATRIX_CSV="${OUTDIR}/dedup_AMR_analytic_matrix.csv"
GENE_LIST="${OUTDIR}/gene_list.mapq${MIN_MAPQ}.txt"

# Threads for samtools
SAMTOOLS_THREADS=8

# Number of samples (must match --array range)
NSAMPLES=60

# Conda env (samtools must be available)
source /data3/ssanchez/miniforge3/etc/profile.d/conda.sh
# conda activate <your_env_with_samtools>

echo "Job started on $(date)"
echo "Running on node(s): $SLURM_NODELIST"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Using BAM dir           : ${AMR_BAM_DIR}"
echo "MAPQ cutoff             : ${MIN_MAPQ}"
echo "Resistome threshold     : ${RESISTOME_THRESHOLD}%"
echo "Output dir              : ${OUTDIR}"

########################
# GET SAMPLE FOR THIS ARRAY TASK
########################

if [[ ! -f "${SAMPLE_LIST_FILE}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

LINE_NO="${SLURM_ARRAY_TASK_ID}"
SAMPLE="$(sed -n "${LINE_NO}p" "${SAMPLE_LIST_FILE}")"

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample name on line ${LINE_NO} in ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

echo "Processing sample: ${SAMPLE}"

########################
# PER-SAMPLE COVERAGE + GENE COUNTS (MAPQ FILTERED + 80% COVERAGE)
########################

BAM="${AMR_BAM_DIR}/${SAMPLE}${BAM_SUFFIX}"

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found for sample ${SAMPLE}: ${BAM}" >&2
    exit 1
fi

# Coverage file (fraction per gene)
COV_FILE="${OUTDIR}/${SAMPLE}.mapq${MIN_MAPQ}.coverage.tsv"

# Counts file
SAMPLE_COUNTS="${OUTDIR}/${SAMPLE}.mapq${MIN_MAPQ}.gene_counts.tsv"
TMP_COUNTS="${SAMPLE_COUNTS}.tmp"

echo "  BAM: ${BAM}"
echo "  Coverage file: ${COV_FILE}"
echo "  Output counts: ${SAMPLE_COUNTS}"

########## 1) Compute coverage fraction per gene (MAPQ-filtered) ##########

# samtools depth:
#  -aa : report all positions, including zero coverage
#  -G 2048 : add SUPPLEMENTARY (2048) to default exclude flags
#            default excluded: UNMAP,SECONDARY,QCFAIL,DUP
#  -Q MIN_MAPQ : min mapping quality
samtools depth -aa -@ "${SAMTOOLS_THREADS}" -G 2048 -Q "${MIN_MAPQ}" "${BAM}" \
    | awk 'BEGIN{OFS="\t"}
           {
               gene = $1; depth = $3;
               len[gene]++;               # total positions (gene length)
               if (depth > 0) cov[gene]++ # positions covered
           }
           END {
               for (g in len) {
                   frac = (cov[g] / len[g]) * 100;
                   printf "%s\t%.2f\n", g, frac;
               }
           }' \
    > "${COV_FILE}"

########## 2) Count reads per gene, keep only genes with coverage ≥ threshold ##########

# Write header to temporary counts file
echo -e "gene_accession\tcount" > "${TMP_COUNTS}"

# samtools view:
#  -F 2308 : exclude unmapped(4), secondary(256), supplementary(2048)
#  -q MIN_MAPQ : enforce mapping quality
# Column 3 ($3) is RNAME = reference / gene accession.
samtools view -@ "${SAMTOOLS_THREADS}" -F 2308 -q "${MIN_MAPQ}" "${BAM}" \
    | awk '{print $3}' \
    | LC_ALL=C sort \
    | uniq -c \
    | awk -v thresh="${RESISTOME_THRESHOLD}" '
        NR==FNR {
            cov[$1] = $2;     # COV_FILE: gene\tcoverage_fraction
            next;
        }
        {
            # From uniq -c: $1 = count, $2 = gene
            count = $1;
            gene  = $2;
            if (cov[gene] + 0 >= thresh)
                printf "%s\t%d\n", gene, count;
        }' "${COV_FILE}" - \
    >> "${TMP_COUNTS}"

# Move atomically to final name
mv "${TMP_COUNTS}" "${SAMPLE_COUNTS}"

echo "Finished counting (with coverage filter) for sample ${SAMPLE} at $(date)"

########################
# MERGE LOGIC (ONLY LAST TASK DOES IT)
########################

# Let the last array task (ID == NSAMPLES) handle the merge
if [[ "${SLURM_ARRAY_TASK_ID}" -ne "${NSAMPLES}" ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID}: counts done, merge will be done by task ${NSAMPLES}."
    exit 0
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: will perform the final merge after all counts are ready."

# Build expected count file list from sample_names.txt
declare -a EXPECTED_COUNTS
while read -r s; do
    [[ -z "${s}" ]] && continue
    [[ "${s}" =~ ^# ]] && continue
    EXPECTED_COUNTS+=( "${OUTDIR}/${s}.mapq${MIN_MAPQ}.gene_counts.tsv" )
done < "${SAMPLE_LIST_FILE}"

# Wait until all expected count files exist
echo "Waiting for all ${#EXPECTED_COUNTS[@]} count files to appear..."
while true; do
    missing=0
    for f in "${EXPECTED_COUNTS[@]}"; do
        if [[ ! -f "${f}" ]]; then
            missing=1
            break
        fi
    done
    if [[ "${missing}" -eq 0 ]]; then
        echo "All count files found at $(date). Proceeding to merge."
        break
    fi
    echo "Not all count files present yet, sleeping 60 seconds..."
    sleep 60
done

########################
# COLLECT COUNT FILES IN ORDER OF sample_names.txt
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
        echo "WARNING: counts file not found for sample ${s}: ${counts} (skipping this sample)" >&2
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

# Concatenate all per-sample counts, skipping header line in each
for f in "${COUNT_FILES[@]}"; do
    if [[ -s "${f}" ]]; then
        tail -n +2 "${f}" >> "${TMP_ALL}"
    fi
done

if [[ ! -s "${TMP_ALL}" ]]; then
    echo "WARNING: No genes found (only headers in all count files). Writing empty matrix."
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
# MERGE INTO WIDE MATRIX: gene_accession x samples
########################

echo "[$(date)] Merging per-sample counts into matrix -> ${MATRIX_CSV}"

# Build space-separated sample name string in the same order as COUNT_FILES
SAMPLE_LIST_STR=""
for s in "${SAMPLES[@]}"; do
    SAMPLE_LIST_STR+="${s} "
done

awk -v samples="${SAMPLE_LIST_STR}" '
BEGIN {
    FS = OFS = "\t";
    split(samples, sname, " ");
    file_idx = 0;
}

# First input file: gene list
NR == FNR {
    genes[++ng] = $1;
    next;
}

# For subsequent files, detect start of each file by FNR == 1 (header line)
FNR == 1 && NR != 1 {
    file_idx++;
    # skip header "gene_accession\tcount"
    next;
}

# Per-sample count files: gene_accession \t count
{
    gene = $1;
    count = $2 + 0;
    acc[gene, file_idx] = count;
}

END {
    sep = ",";
    nfiles = file_idx;

    # Header
    printf "gene_accession";
    for (i = 1; i <= nfiles; i++) {
        if (sname[i] == "") sname[i] = "sample" i;
        printf "%s%s", sep, sname[i];
    }
    printf "\n";

    # Rows
    for (g = 1; g <= ng; g++) {
        gene = genes[g];
        printf "%s", gene;
        for (i = 1; i <= nfiles; i++) {
            key = gene SUBSEP i;
            val = (key in acc ? acc[key] : 0);
            printf "%s%s", sep, val;
        }
        printf "\n";
    }
}
' "${GENE_LIST}" "${COUNT_FILES[@]}" > "${MATRIX_CSV}"

echo "[$(date)] Finished merge."
echo "Final matrix (CSV): ${MATRIX_CSV}"
echo "Dimensions: ${N_GENES} genes x ${#SAMPLES[@]} samples (plus header)"

