#!/usr/bin/env python3
# Verbatim Python block from provenance/original_scripts/bracken.sh
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
