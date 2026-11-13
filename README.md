# Coastal One Health AMR surveillance

A complete Nextflow DSL2 workflow that downloads the PRJEB75359 paired raw reads
from ENA and generates resistance measurements, taxonomic affiliations,
statistical results, supplementary tables and figures from the paper's methods.
Figure 1A is excluded. Historical scripts and historical numerical results are
not needed by the workflow.

The bundled input files in `data/` select 60 filter libraries from 20 biological site–date
samples. The three filters are pooled before abundance normalisation; they are
not treated as independent biological replicates.

## Run

Install Nextflow >=25.04.8, Java 17–21, Python >=3.10, Git and Conda first.
The workflow uses **Conda environments**, with no Docker/Singularity container
configuration. With `-profile conda`, Nextflow automatically creates and caches
the environments declared in `envs/reads.yml` and `envs/analysis.yml`, installing
the analysis tools and R packages on first use. You do not need to install BWA,
Bowtie2, SAMtools, Kraken2, Bracken, fastp, FastQC, MultiQC or R separately.
Internet access and sufficient space for Conda environments are required.
Without the Conda profile, the tools must already be available on PATH.

Keep `data/carpdm_genes.fna.gz` in the repository. It is the default CARPDM
reference and is decompressed automatically before indexing. Run from this
repository:

```bash
nextflow run main.nf -profile conda \
  --kraken_db /absolute/path/k2_standard_20251015 \
  --taxonomy_dir /absolute/path/matching_taxdump \
  --outdir results -resume
```

For Slurm use `-profile conda,slurm`; set your partition/account in a local
Nextflow configuration. Kraken2 requests 128 GB by default; adjust with
`--kraken_memory '64 GB'` or a suitable value for your database/host.
The Slurm profile keeps reference preparation and network downloads on the
controller host, while computational stages run as cluster jobs.

Required biological reference inputs:

* The bundled **CARPDM nucleotide FASTA** at `data/carpdm_genes.fna.gz`, with
  unique reference identifiers. Override with `--carpdm_fasta /path/to/file.fna.gz`
  if needed; plain FASTA also works. Gene labels are read from the final
  pipe-delimited field; editable class annotations are in
  `data/carpdm_gene_classes.tsv`. Unmapped labels remain explicit.
* A **Kraken2 database** and its Bracken `database150mers.kmer_distrib` file.
* The corresponding frozen **NCBI nodes.dmp and names.dmp**. The workflow checks
  the database root or its `taxonomy/` subdirectory automatically; supply
  `--taxonomy_dir` when those files are elsewhere.

MEGARes 3.00, crAssphage NC_024711.1 and the pinned AMR++ SNP verification
resources are retrieved automatically. To reuse local resources, set
`--megares_fasta`, `--crass_fasta` and/or `--snp_bundle`. The latter must contain
`SNP_Verification.py`, `config.ini` and `data/SNPinfo.fasta` plus its Python modules.
Resource checksums and runtime tool versions are recorded. Use the same frozen
Kraken2/taxonomy resource set across runs.

## ENA identifiers and canonical sample IDs

`data/filters.tsv` is the explicit crosswalk: `run_accession` identifies the
ENA run and `filter_id` is your canonical `site_replicate_date` identifier.
For example, `ERR13120681` maps to `A_1_290321`. The resolver retrieves the ENA
FASTQ URLs and MD5 checksums using the run accession, then passes the canonical
filter ID into all downstream processes. Downloaded files initially keep the
ENA run name; the merge stage exposes them as `A_1_290321_R1.fastq.gz` and
`A_1_290321_R2.fastq.gz`, even for a single run. Outputs carry `filter_id`.
The `biosample` column is the internal site/date pooling key (for example,
`A_2021-03-29`), not an ENA BioSample accession. Three filters are pooled into
one such biological sample. Multiple runs per filter can be separated by `;`.
The live run-to-filter download crosswalk is saved in
`results/provenance/ena/downloads.tsv`. The bundled ENA snapshot documents the
mapping; it does not replace live retrieval of URLs and checksums.

## Analysis stages

1. Resolve the explicit ENA run crosswalk, require paired FASTQs and validate
   downloaded MD5 checksums. Single-run filters use symlinks; multiple runs are
   concatenated in matched mate order.
2. FastQC, fastp (Q20, at most 20% unqualified bases, at most five Ns, length
   >=50 bp) and MultiQC.
3. Broad profiling: BWA-MEM alignment to MEGARes, SAMtools duplicate removal,
   MAPQ >=30, exact reference breadth >=80%, and pinned AMR++ SNP verification.
   These run as native Nextflow stages. No nested Nextflow controller is needed.
4. CARPDM Bowtie2 screening: very-sensitive end-to-end alignment, proper pairs,
   MAPQ >=30, breadth >=80% and full-reference mean depth >=5. Relaxed screening
   removes only the depth threshold. Overlapping read ends count separately.
5. Paired Kraken2 classification, the paper's custom 0.2 k-mer support rule,
   frozen rank-aware annotation and Bracken species re-estimation (150 bp,
   minimum ten reads). Unique query/reference associations provide exploratory
   taxonomic affiliations and completeness summaries.
6. crAssphage mapping/coverage, biological-sample pooling, RPM and reference-level
   RPKM, recurrence and filter-support summaries.
7. Hellinger/Bray–Curtis ordinations, date-blocked PERMANOVA and dispersion,
   Friedman and exact paired tests, sequencing-depth correlations, sampling/
   length sensitivities, and partial RDAs conditioned on site and date.
8. Batch PDF/PNG figures and TSV tables. Persistent above-median affiliation
   networks are drawn directly in R; Cytoscape node/edge import tables are also
   exported. An open Cytoscape application is unnecessary.

## Results

`results/manuscript/` contains:

* `figures/`: Figure 1B, Figure 2 and Figures S1–S8. Figure S4 has separate gene
  heatmap and shared-date overlap PDFs.
* `tables/`: Table S1 and S4–S12 components, additional community profiles, network
  nodes/edges and diversity outputs. TSV components preserve different schemas
  without forcing them into one oversized workbook.
* `plot_data/`: plotted values, permutation matrices and environmental models.
* `provenance/`, `sessionInfo.txt` and `INTERPRETATION.txt`.

QC reports and Nextflow execution/trace/DAG reports are written separately under
`results/`. Raw/QC reads and BAMs are retained in `work/` for resumption. Keep
that directory while using `-resume`. Plan storage from the ENA `fastq_bytes`
metadata plus the additional QC/BAM/intermediate files.

Tables S2 and S3 are not generated. Software versions, R packages and reference
checksums are retained as provenance. Environmental measurements are supplied
in `data/filters.tsv`; they cannot be inferred from sequencing reads. Results may differ from
the paper because of database/tool revisions, exact coverage screening, corrected
processing and permutation sampling. No historical result is inserted into a
computed output. See [METHODS.md](docs/METHODS.md) and
[OUTPUTS.md](docs/OUTPUTS.md) for definitions and interpretation.

## Validation

The workflow wiring, a synthetic paired-read tool chain, a known SNP/wild-type
mixture, pooling/taxonomy tests and a complete synthetic 60-filter R report are
covered by tests. The full multi-terabyte ENA study has not been executed here.

```bash
# With the reads environment active:
python -m unittest discover -s tests -v
bash tests/read_smoke.sh /absolute/path/scratch/read-smoke

# Wiring only; uses stubs and produces no scientific results:
nextflow run main.nf -profile test -stub-run --outdir results/stub
```

The test suite's known-allele check uses `AMR_SNP_BUNDLE=/path/to/pinned/verifier`.
See [VALIDATION.md](docs/VALIDATION.md) for tested scope and limits.

