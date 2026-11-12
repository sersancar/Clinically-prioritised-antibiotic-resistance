#!/bin/bash
#SBATCH --job-name=Bracken_taxonomy
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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

# Directory where *filtered* Kraken outputs live (from kraken2Filter.sh)
# Files are assumed to be: <sample>.kraken.score0.2.out
KRAKEN_FILTER_DIR="/data3/ssanchez/Kraken_Taxonomy"
KRAKEN_SCORE=0.2

# Kraken2 / Bracken database directory
KRAKEN_DB="/data3/ssanchez/kraken2DB"
BRACKEN_DB="${KRAKEN_DB}"

# Output directory for Bracken outputs
OUTDIR="/data3/ssanchez/Bracken_taxonomy"
mkdir -p "${OUTDIR}"

# Bracken settings
BRACKEN_LEVEL="S"        # S=species, G=genus, etc.
BRACKEN_MIN_READS=10     # -t
READ_LEN=150             # must match your read length / kmer_distrib

# NCBI taxdump directory for taxonkit (your DB has nodes.dmp/names.dmp in root)
export TAXONKIT_DB="${KRAKEN_DB}"

# Conda env (python, bracken, taxonkit must be available)
source /data3/ssanchez/miniforge3/etc/profile.d/conda.sh
# conda activate <your_env_with_bracken_taxonkit>

command -v python   >/dev/null 2>&1 || { echo "ERROR: python not in PATH" >&2; exit 1; }
command -v bracken  >/dev/null 2>&1 || { echo "ERROR: bracken not in PATH" >&2; exit 1; }
command -v taxonkit >/dev/null 2>&1 || { echo "ERROR: taxonkit not in PATH" >&2; exit 1; }

# Wait until main directory is mounted
while [ ! -d "/data3/ssanchez" ]; do
    echo "Directory /data3/ssanchez does not exist. Waiting..."
    sleep 5
done

echo "Job started on $(date)"
echo "Running on node(s): $SLURM_NODELIST"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Filtered Kraken dir : ${KRAKEN_FILTER_DIR}"
echo "Kraken/Bracken DB   : ${KRAKEN_DB}"
echo "Output dir          : ${OUTDIR}"
echo "Bracken level       : ${BRACKEN_LEVEL}"
echo "Kraken score cutoff : ${KRAKEN_SCORE}"
echo "Read length         : ${READ_LEN}"
echo "TAXONKIT_DB         : ${TAXONKIT_DB}"

########################
# GET SAMPLE NAME FROM ARRAY
########################

if [[ ! -f "${SAMPLE_LIST_FILE}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

LINE_NO="${SLURM_ARRAY_TASK_ID}"
SAMPLE="$(sed -n "${LINE_NO}p" "${SAMPLE_LIST_FILE}" | tr -d '\r' | xargs)"

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample name on line ${LINE_NO} in ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

echo "Processing sample: ${SAMPLE}"

KRAKEN_FILTERED="${KRAKEN_FILTER_DIR}/${SAMPLE}.kraken.score${KRAKEN_SCORE}.out"
if [[ ! -f "${KRAKEN_FILTERED}" ]]; then
    echo "ERROR: Filtered Kraken file not found: ${KRAKEN_FILTERED}" >&2
    exit 1
fi

PREFIX="${OUTDIR}/${SAMPLE}"

BRACKEN_INPUT_REPORT="${PREFIX}.bracken_input.score${KRAKEN_SCORE}.report"
RAW_BRACKEN="${PREFIX}.bracken.${BRACKEN_LEVEL}.tsv"
FINAL_TAX="${PREFIX}.bracken.${BRACKEN_LEVEL}.taxonomy.tsv"

########################
# 1) BUILD A PROPER KRAKEN REPORT FROM FILTERED OUTPUT
########################

echo "[$(date)] Building Kraken-style report from filtered Kraken output -> ${BRACKEN_INPUT_REPORT}"

python - "${KRAKEN_FILTERED}" "${KRAKEN_DB}" "${BRACKEN_INPUT_REPORT}" << 'EOF'
import sys, os, re, collections

if len(sys.argv) != 4:
    sys.stderr.write("Usage: python build_report.py <kraken_filtered> <db_dir> <report_out>\n")
    sys.exit(1)

kraken_file, db_dir, report_path = sys.argv[1:4]

taxid_re = re.compile(r"\(taxid\s+(\d+)\)")

# 1) Count reads directly assigned to each taxid
direct = collections.Counter()
total_reads = 0

with open(kraken_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if not parts:
            continue
        if parts[0] != "C":
            continue
        if len(parts) < 3:
            continue
        name_field = parts[2]
        m = taxid_re.search(name_field)
        if not m:
            continue
        taxid = int(m.group(1))
        if taxid <= 0:
            continue
        direct[taxid] += 1
        total_reads += 1

if total_reads == 0:
    open(report_path, "w").close()
    sys.exit(0)

# 2) Load NCBI taxonomy from Kraken DB (nodes.dmp)
nodes_path = os.path.join(db_dir, "nodes.dmp")
if not os.path.exists(nodes_path):
    sys.stderr.write(f"ERROR: nodes.dmp not found in {db_dir}\n")
    sys.exit(1)

parents = {}
ranks = {}

with open(nodes_path) as f:
    for line in f:
        fields = [x.strip() for x in line.split("|")]
        if len(fields) < 3:
            continue
        try:
            taxid = int(fields[0])
            parent = int(fields[1])
        except ValueError:
            continue
        rank = fields[2]
        parents[taxid] = parent
        ranks[taxid] = rank

# 3) Compute clade counts by propagating direct counts up the tree
clade = collections.Counter()
active = set()

for taxid, cnt in direct.items():
    cur = taxid
    visited = set()
    while True:
        if cur in visited:
            break
        visited.add(cur)
        clade[cur] += cnt
        active.add(cur)
        parent = parents.get(cur)
        if parent is None or parent == cur:
            break
        cur = parent

# Ensure root (taxid 1) is present if available
if 1 in parents:
    active.add(1)
    if clade.get(1, 0) == 0:
        clade[1] = total_reads

# 4) Build children adjacency for the active sub-tree
children = collections.defaultdict(list)
for child in active:
    parent = parents.get(child)
    if parent is None or parent == child:
        continue
    if parent in active:
        children[parent].append(child)

# 5) Load scientific names for active taxa
names = {}
names_path = os.path.join(db_dir, "names.dmp")
if not os.path.exists(names_path):
    sys.stderr.write(f"ERROR: names.dmp not found in {db_dir}\n")
    sys.exit(1)

with open(names_path) as f:
    for line in f:
        fields = [x.strip() for x in line.split("|")]
        if len(fields) < 4:
            continue
        try:
            taxid = int(fields[0])
        except ValueError:
            continue
        if taxid not in active:
            continue
        name_txt = fields[1]
        name_class = fields[3]
        if name_class == "scientific name":
            names[taxid] = name_txt

def rank_code(taxid: int) -> str:
    if taxid == 1:
        return "R"
    r = ranks.get(taxid, "")
    if r == "superkingdom":
        return "D"
    if r == "kingdom":
        return "K"
    if r in ("phylum", "division"):
        return "P"
    if r == "class":
        return "C"
    if r == "order":
        return "O"
    if r == "family":
        return "F"
    if r == "genus":
        return "G"
    if r.startswith("species"):
        return "S"
    return "-"

def sort_children(taxids):
    return sorted(taxids, key=lambda t: (-clade[t], t))

root_taxid = 1 if 1 in active else max(active, key=lambda t: clade[t])

with open(report_path, "w") as out:
    def emit(taxid: int, depth: int):
        clade_cnt = clade.get(taxid, 0)
        direct_cnt = direct.get(taxid, 0)
        pct = (clade_cnt * 100.0 / total_reads) if total_reads > 0 else 0.0
        code = rank_code(taxid)
        name = names.get(taxid, f"taxid_{taxid}")
        indent = " " * (depth * 2)
        out.write(f"{pct:.2f}\t{clade_cnt}\t{direct_cnt}\t{code}\t{taxid}\t{indent}{name}\n")
        for child in sort_children(children.get(taxid, [])):
            emit(child, depth + 1)

    emit(root_taxid, 0)
EOF

# If report is empty (no reads), do not run Bracken
if [[ ! -s "${BRACKEN_INPUT_REPORT}" ]]; then
    echo "No classified reads after Kraken score filter for sample ${SAMPLE}. Writing empty Bracken outputs."
    echo -e "name\ttaxid\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads" > "${RAW_BRACKEN}"
    echo -e "taxid\tname\tlevel\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TAX}"
    exit 0
fi

########################
# 2) RUN BRACKEN
########################

echo "[$(date)] Running Bracken at level ${BRACKEN_LEVEL}, read length ${READ_LEN} -> ${RAW_BRACKEN}"

bracken \
  -d "${BRACKEN_DB}" \
  -i "${BRACKEN_INPUT_REPORT}" \
  -o "${RAW_BRACKEN}" \
  -r "${READ_LEN}" \
  -l "${BRACKEN_LEVEL}" \
  -t "${BRACKEN_MIN_READS}"

N_BRACKEN=$(wc -l < "${RAW_BRACKEN}" || echo 0)
echo "  [INFO] Bracken rows (incl. header): ${N_BRACKEN}"

if [[ "${N_BRACKEN}" -le 1 ]]; then
    echo "WARNING: Bracken output has only header. Writing empty taxonomy table."
    echo -e "taxid\tname\tlevel\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TAX}"
    exit 0
fi

########################
# 3) TAXID -> CANONICAL RANKS (TAXONKIT REFORMAT) AND FINAL TAXONOMY TABLE
########################

TAXIDS_TXT="${PREFIX}.bracken.${BRACKEN_LEVEL}.taxids.txt"
TAXMAP="${PREFIX}.bracken.${BRACKEN_LEVEL}.taxid_ranks.tsv"

echo "[$(date)] Building taxid list from Bracken -> ${TAXIDS_TXT}"
cut -f2 "${RAW_BRACKEN}" | tail -n +2 | sort -u > "${TAXIDS_TXT}"

echo "[$(date)] Building taxid -> canonical ranks (Domain..Species) -> ${TAXMAP}"

# Important findings:
# - taxonkit reformat output does NOT automatically include taxid unless it is part of input line.
# - {k} (kingdom) is often blank in NCBI because many lineages use superkingdom instead.
# Therefore we use {d} (domain/superkingdom) instead of {k}.
taxonkit reformat -I 1 -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" "${TAXIDS_TXT}" > "${TAXMAP}.ranks.tmp"

# Prepend taxid column explicitly so join works correctly
echo -e "taxid\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${TAXMAP}"
paste "${TAXIDS_TXT}" "${TAXMAP}.ranks.tmp" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
  >> "${TAXMAP}"
rm -f "${TAXMAP}.ranks.tmp"

echo "[$(date)] Joining Bracken table with taxonomy -> ${FINAL_TAX}"

# Reorder Bracken so that taxid is first column:
# taxid  name  level  kraken_reads  added_reads  new_est_reads  fraction
BRACKEN_REORDER="${RAW_BRACKEN}.reorder"
tail -n +2 "${RAW_BRACKEN}" \
  | awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1, $3, $4, $5, $6, $7}' \
  > "${BRACKEN_REORDER}"

# Sort both on taxid
LC_ALL=C sort -t $'\t' -k1,1 "${BRACKEN_REORDER}" > "${BRACKEN_REORDER}.sorted"
LC_ALL=C tail -n +2 "${TAXMAP}" | sort -t $'\t' -k1,1 > "${TAXMAP}.sorted"

# Join on taxid; keep all Bracken rows even if taxonomy is missing for some taxids
join -t $'\t' -a 1 -e "" \
  -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7,2.8 \
  "${BRACKEN_REORDER}.sorted" \
  "${TAXMAP}.sorted" \
  > "${FINAL_TAX}.joined"

# Final header (Domain instead of Kingdom, per our findings)
echo -e "taxid\tname\tlevel\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction\trank\tDomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > "${FINAL_TAX}"

awk -F'\t' 'BEGIN{OFS="\t"}
{
    taxid   = $1
    if (taxid == "1") next

    name    = $2
    level   = $3
    kreads  = $4
    add     = $5
    est     = $6
    frac    = $7

    domain  = $8
    phylum  = $9
    classs  = $10
    order   = $11
    family  = $12
    genus   = $13
    species = $14

    rank = ""
    if (species != "")      rank = "species"
    else if (genus != "")   rank = "genus"
    else if (family != "")  rank = "family"
    else if (order != "")   rank = "order"
    else if (classs != "")  rank = "class"
    else if (phylum != "")  rank = "phylum"
    else if (domain != "")  rank = "domain"

    print taxid, name, level, kreads, add, est, frac,
          rank, domain, phylum, classs, order, family, genus, species
}' "${FINAL_TAX}.joined" >> "${FINAL_TAX}"

rm -f "${BRACKEN_REORDER}" "${BRACKEN_REORDER}.sorted" "${TAXMAP}.sorted" "${FINAL_TAX}.joined"

echo "[$(date)] Finished Bracken-based taxonomy for sample: ${SAMPLE}"
echo "Bracken report input : ${BRACKEN_INPUT_REPORT}"
echo "Bracken table (raw)  : ${RAW_BRACKEN}"
echo "Bracken+NCBI taxonomy: ${FINAL_TAX}"

