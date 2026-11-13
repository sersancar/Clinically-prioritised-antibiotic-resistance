# Output mapping

| Output prefix | Content |
|---|---|
| TableS1 | Filter metadata, biological correspondence, runtime denominators, environmental values |
| TableS4 | Complete CARPDM per-filter reference measurements, primary/relaxed flags |
| TableS5 | CARPDM biological counts, support, RPM/RPKM under both screens |
| TableS6 | Site/study threshold sensitivity, with target/reference counting units distinguished |
| TableS7 | CARPDM taxonomic affiliations and filter/site/study completeness |
| TableS8 | Broad biological counts, annotations, RPM/RPKM and shared-date gene overlap |
| TableS9 | Broad resistance-associated taxonomy |
| TableS10 | PERMANOVA, dispersion, Friedman, exact contrasts and depth correlations |
| TableS11 | Partial RDAs, predictor audit and tidal context |
| TableS12 | Marker filter measurements, biological RPM and site summaries |

Tables S2 and S3 are not generated. Original table numbering is retained for
the data-derived tables. Software versions, R packages and reference checksums
are saved under `provenance/`, outside the manuscript tables.

`provenance/kraken_database.json` records the automatic database download URLs,
pinned checksum manifest and verified file hashes. The downloaded database
and taxonomy remain in the cached Nextflow task under `work/`.

Additional TSVs retain community profiles, broad per-filter counts, diversity,
network node/edge selections, CARPDM means/intervals and heatmap pseudocounts.
Different components of one supplementary table remain separate TSVs, preserving
numeric types and clear schemas. Figures are PDF and PNG, with a separate PDF
for the S4 Venn overlap. Source data for plots and permutation/model objects are
saved alongside them. Figure 1A is not generated.

Empty screening results generate labelled empty panels. Ordination requires
nonzero broad profiles; invalid or empty biological profiles stop the report
with an explicit explanation. Unknown CARPDM class annotations remain Unassigned,
with exact target/gene identifiers preserved in the tables.
