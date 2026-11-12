#!/usr/bin/env bash
# Use an environment with the reads.yml tools. No ENA reads are downloaded.
set -euo pipefail
project=$(cd "$(dirname "$0")/.." && pwd)
destination=$(realpath -m "${1:?Supply a scratch output directory}")
mkdir -p "$destination"
export PATH="$project/bin:$PATH"
python "$project/tests/make_reads.py" "$destination"
cd "$destination"
mkdir -p qc
fastqc --threads 2 --outdir qc reads_1.fq.gz reads_2.fq.gz
fastp -i reads_1.fq.gz -I reads_2.fq.gz -o clean_1.fq.gz -O clean_2.fq.gz --qualified_quality_phred 20 \
  --unqualified_percent_limit 20 --n_base_limit 5 --length_required 50 --json sample.fastp.json --html qc.html --thread 2
cp sample.fastp.json qc/
multiqc qc --force --outdir multiqc
bwa index megares.fa
bwa mem -t 2 megares.fa clean_1.fq.gz clean_2.fq.gz | samtools view -u -F 2308 - | samtools sort -n -o names.bam -
samtools fixmate -m names.bam fixed.bam
samtools sort -o coordinate.bam fixed.bam
samtools markdup -r coordinate.bam dedup.bam
samtools view -b -q 30 -F 2308 dedup.bam > broad.bam
samtools index broad.bam
bam_metrics.py --bam broad.bam --fasta megares.fa --filter-id sample --catalogue broad --output broad.metrics.tsv
python - <<'PY'
from common import read_tsv,write_tsv
rows=read_tsv('broad.metrics.tsv')
write_tsv('broad.counts.tsv',[dict(r,count=int(r['mapped_read_ends']) if r['primary']=='True' else 0) for r in rows])
assert all(r['primary']=='True' for r in rows),rows
PY
bowtie2-build carpdm.fa index
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -k 1 -x index -1 clean_1.fq.gz -2 clean_2.fq.gz | \
  samtools view -b -f 2 -q 30 -F 2308 - | samtools sort -o carpdm.bam -
samtools index carpdm.bam
bam_metrics.py --bam carpdm.bam --fasta carpdm.fa --filter-id sample --catalogue carpdm --output carpdm.metrics.tsv
mkdir -p db/taxonomy
cp taxonomy/*.dmp db/taxonomy/
kraken2-build --add-to-library genomes.fa --db db --no-masking
kraken2-build --build --db db --threads 2 --no-masking
kraken2 --db db --paired --gzip-compressed --threads 2 --output kraken.tsv clean_1.fq.gz clean_2.fq.gz
link_reads.py --filter-id sample --broad-bam broad.bam --carpdm-bam carpdm.bam --broad-counts broad.counts.tsv \
  --carpdm-metrics carpdm.metrics.tsv --kraken kraken.tsv --taxonomy taxonomy --out links
bracken-build -d db -t 2 -l 150
bracken -d db -i links/community.report -o species.tsv -r 150 -l S -t 10
python - <<'PY'
from common import read_tsv
metrics=read_tsv('carpdm.metrics.tsv');comp=read_tsv('links/completeness.tsv')
assert len(metrics)==2 and all(r['primary']=='True' for r in metrics)
assert all(int(r['eligible_associations'])==122 for r in comp),comp
assert all(int(r['genus_resolved'])==122 for r in comp),comp
assert len(read_tsv('species.tsv'))==2
print('Paired-read QC, mapping, deduplication, coverage, taxonomy linkage and Bracken smoke test passed.')
PY
