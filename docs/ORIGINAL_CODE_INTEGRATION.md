# Original-code integration — 12 November 2026

All fourteen uploaded scripts are archived unchanged under
`provenance/original_scripts/`, with uploaded filenames and SHA256 checksums in
`manifest.json`. Runnable derivatives live in `bin/`, `analysis/` and `adapters/`.

## Integrated stages

* CARPDM now uses the original end-to-end, very-sensitive mapping with
  `--no-mixed --no-discordant -k 1`, and properly paired alignments (`-f 2`).
  crAssphage retains its separate very-sensitive mapping specification.
* Custom Kraken2 scoring comes from `kraken2Filter.sh`. Its score calculation
  is unchanged. Deduplication moves from an in-memory set to external sort;
  tests compare retained records with the original embedded Python block.
* Report reconstruction is the original Python block from `bracken.sh`.
  The workflow explicitly requests Kraken2 `--use-names`, matching the report
  builder's `name (taxid N)` parser. Bracken settings remain S, 150 bp and 10.
  ARG/taxonomy joins and rank-aware completeness remain separate work.
* `prepare_r_inputs.py` translates validated products into the original R
  input schema. It uses runtime fastp read totals and maps freshwater filter
  IDs to their original `F1FW...` form. It requires already-filtered AMR counts
  and existing CARPDM taxonomy tables; it cannot create missing SNP screening.
* `FROM_TABLES` runs the supplied R analysis, plots and optional offline networks
  from prepared tables. It does not represent precomputed inputs as a raw-read
  analysis.

## Portable R changes

Analysis paths use `AREST_INPUT`, `AREST_OUTPUT` and `AREST_RUN`.
`AREST_CARPDM_UNFILTERED` defaults to FALSE. Seed 20261030 and 9999 permutations
are unchanged. Environmental BH adjustment now covers all three predictors
together, matching the final supplement and the old table assembler's
`Q_BH_all_three`; the source analysis split the family by Analysis label.
Raw test statistics and P values are unaffected by that adjustment change.

All plotting paths use CLI arguments, including the two embedded hard-coded
paths. Heatmap helper definitions/type colours were moved before first use.
Four superseded Venn drafts were removed, retaining the last complete block;
earlier drafts referenced variables absent from a clean session. A later generic
Figure 2A overwrite was removed so its customised standalone version agrees
with the combined figure. These changes require visual comparison before release.

The network script's unconditional final RCy3 layout adjustment is now guarded
by explicit Cytoscape enablement and successful connection. Static R output is
the batch default. The Nextflow network process works on a private analysis copy
and requires a frozen Taxon/Domain map. Live Entrez lookup is disabled by default;
outside the reproducible workflow it can be explicitly enabled with
`AREST_ALLOW_LIVE_TAXONOMY=TRUE`.

## Run the original analysis from tables

```bash
nextflow run main.nf -entry FROM_TABLES -profile conda \
  --analysis_inputs /absolute/path/prepared_inputs \
  --broad_taxonomy /absolute/path/broad_taxonomy \
  --domain_map /absolute/path/taxon_domain_map.tsv \
  --carpdm_unfiltered true --outdir results/reanalysis -resume
```

The network arguments are optional and must be supplied together. Enable
`carpdm_unfiltered` only for verified complete all-reference, all-filter metrics.
Inputs use the original schemas:

| Input | Required content |
|---|---|
| `metadata.tsv` | `Sample_name` and environmental aliases recognised in the R script |
| `sample_manifest.tsv` | `SampleID,Site,Date,BioSample` |
| `total_reads.tsv` | `sample,total_reads`; runtime post-QC read ends |
| `dedup_AMR_analytic_matrix2.csv` | `gene_accession` and 60 library count columns |
| `megares_lengths.tsv` | `gene_accession,length`; exact reference IDs |
| `CARD_genes/ID.ARG_metrics.tsv` | `gene,length,cov_bases,breadth,mean_depth,num_reads`; every filter |
| `CARPDM_taxonomy/ID.CARD_ARG_taxonomy.tsv` | Original schema including `ARG_ID,read_count,genus`; header-only for zero calls |

Network input contains one `ID.ARG_taxonomy.tsv` per filter, with original
`ARG_ID,class,read_count` fields. The domain TSV has `Taxon,Domain` and must
cover all eligible taxonomic labels. Critical R modules are checked even when
the original script records a final COMPLETE status.

## Numbering crosswalk for final packaging

| Original output | Final destination |
|---|---|
| Original schematic, not supplied | Figure 1A |
| `Figure1A_CARPDM_sqrt` | Figure 1B |
| `Figure2_combined` | Figure 2 |
| `FigureS5_temporal_CARPDM` | Figure S1 |
| `FigureS10_gene_filter_support` | Figure S2 |
| `Figure1B_taxonomy_bubbles` | Figure S3 |
| `FigureS2_original_layout_20_samples` and separate Venn | Figure S4 |
| `FigureS6_QC_diversity` | Figure S5 |
| Four `FigureS3_[site]` networks | Figure S6 |
| `FigureS7_design_RPM_RPKM` | Figure S7 |
| `FigureS8_pRDA_site_date` | Figure S8 |

Outputs retain original names pending packaging. The legacy assembler is kept
as `arest_final_supplementary_tables_legacy.R` and is not called automatically.
It expects `S1_Metadata`, `S2_previous_studies` and `S6_ARG_taxonomy` sheets and
produces an earlier S1–S11 arrangement. The current supplied workbook uses
S1–S12 and section headers; a simple sheet rename cannot reproduce it.

## Remaining evidence

* Effective AMR++ configuration/source commit and reference/SNP manifest. The
  uploaded BAM count-rebuild script applies MAPQ and rounded percentage breadth
  but contains no SNP gate. Confirm how SNP screening propagated to its matrix.
* Exact reference FASTAs/annotations and frozen taxonomy/Bracken resources.
  A path/URL/checksum manifest is sufficient initially for large resources.
* SAMtools version and which CARPDM coverage branch actually ran. The originals
  choose `coverage` when available and otherwise `depth`. Compatibility of all
  supplied flags with the table's version 1.8 remains untested.
* Later final S1–S12/completeness code or corresponding numeric intermediates,
  if newer revisions supersede the supplied legacy assembler.
* Frozen network domain map and Figure 1A source. The 57-genus candidate list
  is now present in the supplied plotting script and is no longer missing.

The original TaxonKit linkage scripts are retained pending frozen taxonomy and
validation. They have not been replaced with guessed physical host assignments.

## Additional AMR uploads

The later launcher/CARDGenes/ncahill matrix scripts are now archived as well.
The source archive contains fourteen files. The launcher confirms SNP and
deduplication flags; the two matrix versions differ in duplicate-flag counting.
See AMR_EXECUTION.md for the portable launcher and explicit count-policy
workflow, including the unresolved choice of the final matrix policy.
