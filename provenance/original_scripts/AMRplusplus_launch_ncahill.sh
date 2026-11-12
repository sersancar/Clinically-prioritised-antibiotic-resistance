#!/usr/bin/env bash
#SBATCH --job-name=AMR++
#SBATCH --output=/home/ncahill/sergio/logs/%x_%j.out
#SBATCH --error=/home/ncahill/sergio/logs/%x_%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=50G
#SBATCH --time=500:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sergio.sanchezcarrillo@universityofgalway.ie

# Strict mode (we'll relax -u only during activation)
set -e -o pipefail
set -u

ENV_PREFIX="/data4/ncahill/.conda/envs/AMR++_env"
MINICONDA_PREFIX="/data4/ncahill/miniconda3"
export NXF_CLEANUP=false

# --- 1) Activate the env (disable -u briefly because openjdk's hook isn't nounset-safe)
set +u
if [ -f "$MINICONDA_PREFIX/etc/profile.d/conda.sh" ]; then
  source "$MINICONDA_PREFIX/etc/profile.d/conda.sh"
  conda activate "$ENV_PREFIX" || true
else
  eval "$("$MINICONDA_PREFIX/bin/conda" shell.bash hook)" || true
  conda activate "$ENV_PREFIX" || true
fi
set -u

# --- 2) If activation didn't set CONDA_PREFIX, fallback by path
if [ "${CONDA_PREFIX:-}" != "$ENV_PREFIX" ]; then
  echo "WARN: conda activate didn't set CONDA_PREFIX; applying manual PATH fallback."
  export CONDA_PREFIX="$ENV_PREFIX"
  export PATH="$ENV_PREFIX/bin:$PATH"
fi

# --- 3) Force Java 17 from this env; ignore any old overrides
unset JAVA_CMD NXF_JAVA_CMD
export JAVA_HOME="$CONDA_PREFIX"
export NXF_JAVA_HOME="$CONDA_PREFIX"
export PATH="$JAVA_HOME/bin:$PATH"
hash -r

# --- 4) Sanity logs
echo "Job started on $(date)"
echo "Node(s): ${SLURM_NODELIST:-unknown}"
echo "PWD: $(pwd)"
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "which java: $(command -v java)"
java -version

NEXTFLOW_BIN="$CONDA_PREFIX/bin/nextflow"
if [ ! -x "$NEXTFLOW_BIN" ]; then
  echo "ERROR: $NEXTFLOW_BIN not found or not executable."
  echo "Hint: install with: mamba install -n AMR++_env -c conda-forge -c bioconda nextflow openjdk=17"
  exit 1
fi
echo "which nextflow: $NEXTFLOW_BIN"
"$NEXTFLOW_BIN" -version

# --- 5) Paths/inputs
INPUTDIR="/data4/ncahill/fastpOut/batch2"
OUTPUTDIR="/data4/ncahill/AMR++Out/batch2"
AMRDIR="/data4/ncahill/AMRplusplus"
mkdir -p "$OUTPUTDIR"
cd "$AMRDIR"

# --- 6) Run Nextflow from the env explicitly
"$NEXTFLOW_BIN" run main_AMR++.nf \
  -profile local_slurm \
  -resume \
  --snp Y \
  --deduped Y \
  --threads "${SLURM_CPUS_PER_TASK:-36}" \
  --pipeline resistome \
  --reads "${INPUTDIR}/*_clean_{1,2}.fq.gz" \
  --output "${OUTPUTDIR}"

echo "Job finished on $(date)"

# RUN: sbatch AMR++.sh
