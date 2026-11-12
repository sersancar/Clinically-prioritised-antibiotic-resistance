# Clinically prioritised antibiotic resistance genes identify One Health surveillance hotspots along a coastal wastewater-to-beach gradient

Nextflow DSL2 workflow for the One Health Short Communication using ENA study
PRJEB75359. The design comprises 60 filter libraries from 20 site-date samples
across recreational beaches A/C, untreated wastewater discharge B and
agricultural freshwater input D.

**Status: development; not yet a complete reproducibility release.** ENA
retrieval, QC, CARPDM/crAssphage screening and read-weighted pooling are
implemented. The supplied custom Kraken2/Bracken code and R scripts now have
portable execution paths. Effective AMR++ configuration/resources, ARG/taxonomy joins, final
packaging and numerical validation are still required to reproduce all Figures
and Tables S1–S12. Full mode fails early when required adapters are missing.

## Requirements

Nextflow >=25.04.8, Java 17 or later (Java 21 recommended), Python >=3.10 and
Conda for `-profile conda`. See Table S3-derived
`metadata/software_provenance.tsv`. Tool versions are declared in
`envs/reads.yml`; platform-specific dependency resolution is not yet tested.
R dependencies are declared separately in `envs/analysis.yml` and remain
unresolved. See [original-code integration](docs/ORIGINAL_CODE_INTEGRATION.md)
for the archived scripts, portable changes and numbering crosswalk.

The supplied AMR++ launcher confirms SNP screening and deduplication switches.
See [AMR_EXECUTION.md](docs/AMR_EXECUTION.md) for its revision-checked external
launcher and the `FROM_AMR_BAMS` entry point. The two uploaded BAM-count scripts
use different duplicate-flag policies; an explicit policy is required.

## Run

All 60 ENA run mappings are resolved from descriptive sample titles (see
`metadata/ENA_MAPPING.md`). Supply the remaining verified inputs listed in
[INPUTS_NEEDED.md](docs/INPUTS_NEEDED.md), including original references with
SHA256 values in a copy of `assets/reference_lock.example.json`.

```bash
# Validate metadata without downloading reads.
python -c "import sys; sys.path.insert(0,'bin'); from ena_manifest import *; validate_filters(read_tsv('metadata/filters.tsv'), require_runs=False)"

# Resolve accession selection; a verified mapping is mandatory.
python bin/ena_manifest.py --filters metadata/filters.tsv --outdir ena

# Test workflow wiring with synthetic reads; produces no manuscript results.
nextflow run main.nf -profile test --mode screening -stub-run --outdir results/stub

# Real screening, when reference provenance and original options are confirmed.
nextflow run main.nf -profile conda --mode screening \
  --carpdm_fasta /absolute/path/carpdm.fa \
  --crass_fasta /absolute/path/crassphage_NC_024711.1.fasta \
  --reference_lock /absolute/path/reference_lock.json \
  --outdir results/screening -resume

# Full workflow: requires implemented original-code adapters.
nextflow run main.nf -profile conda,slurm --mode full \
  --carpdm_fasta /absolute/path/carpdm.fa \
  --crass_fasta /absolute/path/crassphage_NC_024711.1.fasta \
  --reference_lock /absolute/path/reference_lock.json \
  --kraken_db /absolute/path/k2_standard_20251015 \
  --amr_adapter /absolute/path/amr_original.sh \
  --taxonomy_adapter /absolute/path/taxonomy_original.sh \
  --report_adapter /absolute/path/manuscript_original.sh \
  --outdir results/full -resume
```

Raw data are large. No production reads have been downloaded or analyses
executed as part of preparing this repository. Keep the Nextflow work directory
for resumption; BAMs and QC reads are retained there and are not published by
default. `--mode download` resolves and downloads raw reads into the work
directory only. Execution reports and trace files go under `outdir/provenance`.

## Outputs and methods

The `FROM_TABLES` entry point reruns the original R analysis and plots, plus
optional offline networks, from prepared inputs. See
[ORIGINAL_CODE_INTEGRATION.md](docs/ORIGINAL_CODE_INTEGRATION.md) for its command
and schemas. Outputs retain the source scripts' numbering pending packaging.

Screening produces per-filter all-reference measurements and pooled
`carpdm.biological.tsv` / `crass.biological.tsv`, with separate primary and
relaxed rules. The report adapter must create Figures 1–2 and S1–S8 plus
Tables S1–S12. Table S2 is curated historical evidence; measured environmental
metadata and the schematic are supplied assets. Their provenance is retained
instead of pretending they can be computed from reads.

See [ADAPTERS.md](docs/ADAPTERS.md) for original-code interfaces,
[REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for analytical rules and release
checks, and `assets/expected_outputs.json` for the report contract.

## Development checks

```bash
python -m unittest discover -s tests -v
python -m compileall -q bin tests
nextflow run main.nf -profile test --mode screening -stub-run --outdir results/stub
```

No software licence has yet been selected. Publication requires the authors'
licence choice, authenticated GitHub access, original-code integration and
numerical/graphical validation. This development repository must not yet be
described as the code used to generate the manuscript results.
