# Original-code integration contracts

Adapters are Bash entry points supplied as absolute paths. They must exit
nonzero on failure. Their dependencies must exist in the target environment;
the core conda environment does not provide the original R/Cytoscape stack.
If an adapter calls auxiliary files, use installed absolute paths or package
them in a container: only the adapter file itself is staged by Nextflow.

## AMR++ adapter

```
bash adapter.sh FILTER_ID QC_R1 QC_R2 OUTPUT_DIRECTORY
```

Run the original AMR++ version with duplicate removal and SNP screening enabled,
MAPQ >=30 and breadth >=0.8 per filter. Do not replace SNP screening by plain
Bowtie2 mapping. Required outputs:

* `reference_counts.tsv`: columns `gene_accession,count`, containing original
  exact-reference mapped-read-end counts. Retain any additional annotations.
* `deduplicated.bam`: BAM used for broad-resistome taxonomic associations.
* `provenance.json`: commit, resource checksums, command, tool versions and
  count definition, including all deduplication/SNP settings.

The report adapter must consume the approved schema of the original counts.
No invented converter schema is imposed before reviewing that script.

## Taxonomy adapter

```
bash adapter.sh FILTER_ID RAW_KRAKEN_OUTPUT KRAKEN_DATABASE OUTPUT_DIRECTORY
```

Apply the published custom score rule, retain uncalculable scores, reconstruct
reports and run Bracken at species rank, length 150, threshold 10. Required:
`retained_assignments.tsv`, `bracken_species.tsv`, `provenance.json`. Record
read/query name handling and use the original taxonomy dump for ancestor ranks.

The default adapter now implements scoring/report reconstruction from the
uploaded code. It produces retained assignments and Bracken species estimates;
ARG-associated rank annotations and completeness tables are not yet integrated.

## Manuscript report adapter

```
bash adapter.sh FILTERS AMR_DIR TAXONOMY_DIR SCREEN_BAMS_DIR SCREEN_METRICS_DIR BIOLOGICAL_DIR QC_DIR OUTPUT_DIR
```

Directory inputs contain all filters; screen BAMs are named
`FILTER.carpdm.bam` / `FILTER.crass.bam`. Metrics preserve all references,
including zero detections. Biological TSVs have both primary and relaxed
screening results. The QC directory contains the fastp JSON for every filter.

Integrate the final original analysis and plotting scripts, including their
taxonomy joins. Export every file specified in `assets/expected_outputs.json`.
Table S1/S2/S3 include measured/curated inputs; pass their verified source files
to the adapter via its installed resource directory. The output contract
checks nonempty files; numerical and graphical agreement require regression
validation against the supplied workbook and figures.

For full production runs, package each adapter's environment and record its
image digest or explicit lockfile. Do not launch a second unpinned Nextflow
workflow invisibly inside an adapter; record its revision and command if the
original AMR++ workflow must run as an external, separately resumed stage.
