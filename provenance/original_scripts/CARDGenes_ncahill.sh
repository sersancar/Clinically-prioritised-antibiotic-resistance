#!/bin/bash
#SBATCH --job-name=CARDGenes
#SBATCH --output=/home/ncahill/sergio/logs/%x_%A_%a.out
#SBATCH --error=/home/ncahill/sergio/logs/%x_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-68
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sergio.sanchezcarrillo@universityofgalway.ie

set -euo pipefail

# ----------------------------
# Conda environment (edit)
# ----------------------------
source /data4/ncahill/miniconda3/etc/profile.d/conda.sh
conda activate samtools_new

# ----------------------------
# Parameters for "presence"
# ----------------------------
MIN_MAPQ=30
MIN_BREADTH=0.80      # fraction of gene length covered by >=1 read
MIN_MEAN_DEPTH=5      # mean depth across full gene length

# Bowtie2 mapping behavior for strict presence evidence
BOWTIE2_MODE=(--end-to-end --very-sensitive --no-mixed --no-discordant -k 1)

# ----------------------------
# Paths
# ----------------------------
INPUTDIR="/data4/ncahill/fastpOut/batch2"
OUTPUTDIR="/data4/ncahill/card_genes"

# Existing Bowtie2 index prefix (as per your ls: card_genes_bt2.1.bt2 etc.)
INDEX_PREFIX="${OUTPUTDIR}/card_genes_bt2"

# Reference FASTA (nucleotide)
REFERENCEFILE="${OUTPUTDIR}/card_genes.fna"

SAMPLE_LIST_FILE="/home/ncahill/sergio/samplesBatch2.txt"

ALIGNDIR="${OUTPUTDIR}/alignments"
STATSDIR="${OUTPUTDIR}/stats"
mkdir -p "${ALIGNDIR}" "${STATSDIR}"

THREADS="${SLURM_CPUS_PER_TASK}"

# ----------------------------
# Debug info
# ----------------------------
echo "Job started on $(date)"
echo "Job: ${SLURM_JOB_ID}  ArrayTask: ${SLURM_ARRAY_TASK_ID}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Workdir: $(pwd)"
echo "CPUs: ${THREADS}"

# ----------------------------
# Basic tool checks
# ----------------------------
command -v bowtie2 >/dev/null 2>&1 || { echo "ERROR: bowtie2 not in PATH" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not in PATH" >&2; exit 1; }

# ----------------------------
# Get sample name
# ----------------------------
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

# ----------------------------
# Fixed read naming (no pattern checks)
# ----------------------------
READ1="${INPUTDIR}/${SAMPLE}_clean_1.fq.gz"
READ2="${INPUTDIR}/${SAMPLE}_clean_2.fq.gz"

if [[ ! -f "${READ1}" ]]; then echo "ERROR: Missing ${READ1}" >&2; exit 1; fi
if [[ ! -f "${READ2}" ]]; then echo "ERROR: Missing ${READ2}" >&2; exit 1; fi

# ----------------------------
# Verify index exists (no rebuilding)
# ----------------------------
if [[ ! -f "${INDEX_PREFIX}.1.bt2" && ! -f "${INDEX_PREFIX}.1.bt2l" ]]; then
  echo "ERROR: Bowtie2 index not found. Expected ${INDEX_PREFIX}.1.bt2 (or .bt2l)" >&2
  exit 1
fi

# Reference FASTA is used only for fallback coverage calculation and optional faidx.
if [[ ! -f "${REFERENCEFILE}" ]]; then
  echo "ERROR: Reference FASTA not found: ${REFERENCEFILE}" >&2
  exit 1
fi

# Create FASTA index (lightweight, helps fallback + transparency)
if [[ ! -f "${REFERENCEFILE}.fai" ]]; then
  samtools faidx "${REFERENCEFILE}"
fi

# ----------------------------
# Outputs
# ----------------------------
BOWTIE2_LOG="${STATSDIR}/${SAMPLE}.bowtie2.log"

BAM_SORTED="${ALIGNDIR}/${SAMPLE}.sorted.bam"
BAM_FILT_SORTED="${ALIGNDIR}/${SAMPLE}.filt.sorted.bam"

FLAGSTAT_TXT="${STATSDIR}/${SAMPLE}.flagstat.txt"
IDXSTATS_TXT="${STATSDIR}/${SAMPLE}.idxstats.txt"

COV_RAW="${STATSDIR}/${SAMPLE}.samtools_coverage.tsv"
METRICS_TSV="${STATSDIR}/${SAMPLE}.ARG_metrics.tsv"
PRESENT_TSV="${STATSDIR}/${SAMPLE}.ARG_present.tsv"

# ----------------------------
# 1) Align -> sorted BAM
# ----------------------------
echo "Aligning with bowtie2..."
bowtie2 \
  "${BOWTIE2_MODE[@]}" \
  -x "${INDEX_PREFIX}" \
  -1 "${READ1}" -2 "${READ2}" \
  -p "${THREADS}" \
  2> "${BOWTIE2_LOG}" \
| samtools view -@ "${THREADS}" -b - \
| samtools sort -@ "${THREADS}" -o "${BAM_SORTED}" -

samtools index -@ "${THREADS}" "${BAM_SORTED}"

# ----------------------------
# 2) Filter for strong evidence (presence-focused)
#    - properly paired (-f 2)
#    - MAPQ >= MIN_MAPQ
#    - remove secondary (256) and supplementary (2048)
# ----------------------------
echo "Filtering BAM for high-confidence evidence (MAPQ>=${MIN_MAPQ}, proper pairs)..."
samtools view -@ "${THREADS}" -b \
  -f 2 -q "${MIN_MAPQ}" \
  -F 256 -F 2048 \
  "${BAM_SORTED}" \
| samtools sort -@ "${THREADS}" -o "${BAM_FILT_SORTED}" -

samtools index -@ "${THREADS}" "${BAM_FILT_SORTED}"

# ----------------------------
# 3) Basic mapping stats
# ----------------------------
samtools flagstat -@ "${THREADS}" "${BAM_FILT_SORTED}" > "${FLAGSTAT_TXT}"
samtools idxstats "${BAM_FILT_SORTED}" > "${IDXSTATS_TXT}"

# ----------------------------
# 4) Per-gene breadth/depth metrics + presence calls
#    Prefer: samtools coverage
#    Fallback: samtools depth + reference .fai
# ----------------------------
echo -e "gene\tlength\tcov_bases\tbreadth\tmean_depth\tnum_reads" > "${METRICS_TSV}"
echo -e "gene\tlength\tcov_bases\tbreadth\tmean_depth\tnum_reads" > "${PRESENT_TSV}"

if samtools coverage -h >/dev/null 2>&1; then
  echo "Computing gene coverage with: samtools coverage"
  samtools coverage "${BAM_FILT_SORTED}" > "${COV_RAW}"

  # samtools coverage columns (typical):
  # rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
  awk -v OFS="\t" -v mb="${MIN_BREADTH}" -v md="${MIN_MEAN_DEPTH}" '
    BEGIN {}
    $0 ~ /^#/ {next}   # skip header lines like #rname ...
    NF < 7 {next}
    {
      gene=$1; start=$2; end=$3; numreads=$4; covbases=$5; meandepth=$7;
      len=end-start+1;
      breadth=(len>0 ? covbases/len : 0);
      print gene,len,covbases,breadth,meandepth,numreads
      if (breadth>=mb && meandepth>=md) {
        print gene,len,covbases,breadth,meandepth,numreads >> "'"${PRESENT_TSV}"'"
      }
    }
  ' "${COV_RAW}" >> "${METRICS_TSV}"

else
  echo "samtools coverage not available; falling back to samtools depth (-aa) + FASTA .fai"

  # Depth output: rname pos depth
  # We compute:
  #  - covered positions (depth>0)
  #  - sum depth across all positions
  #  - mean depth across full gene length
  #  - breadth = covered/length
  #
  # Note: -aa prints all positions for all references in BAM header, which is fine here
  samtools depth -aa -d 0 "${BAM_FILT_SORTED}" \
  | awk -v OFS="\t" -v fai="${REFERENCEFILE}.fai" -v mb="${MIN_BREADTH}" -v md="${MIN_MEAN_DEPTH}" '
      BEGIN {
        # load lengths from .fai
        while ((getline < fai) > 0) {
          len[$1]=$2
        }
        close(fai)
      }
      {
        g=$1; d=$3
        sumd[g]+=d
        if (d>0) covpos[g]++
      }
      END {
        for (g in len) {
          L=len[g]
          cb=(covpos[g]+0)
          breadth=(L>0 ? cb/L : 0)
          meandepth=(L>0 ? sumd[g]/L : 0)
          # num_reads not directly available from depth; set NA
          print g,L,cb,breadth,meandepth,"NA"
          if (breadth>=mb && meandepth>=md) {
            print g,L,cb,breadth,meandepth,"NA" >> "'"${PRESENT_TSV}"'"
          }
        }
      }
    ' >> "${METRICS_TSV}"
fi

# ----------------------------
# Done
# ----------------------------
echo "Done: ${SAMPLE}"
echo "Filtered BAM: ${BAM_FILT_SORTED}"
echo "Metrics:      ${METRICS_TSV}"
echo "Present ARGs: ${PRESENT_TSV}"
echo "Job finished on $(date)"

