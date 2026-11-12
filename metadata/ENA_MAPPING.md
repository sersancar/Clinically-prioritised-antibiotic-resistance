# Verified ENA mapping

On 12 November 2026 the ENA `read_run` file report for PRJEB75359 returned
60 paired runs. The unchanged TSV, including descriptive sample titles and MD5
checksums, is saved as `ena_PRJEB75359_2026-11-12.tsv`. It contains no coordinates.

Each title maps uniquely to one canonical manuscript filter ID. Examples:

| ENA sample title | Manuscript filter | Run |
|---|---|---|
| B1Seawater290321 | B_1_290321 | ERR13120693 |
| 1Freshwater140421 | D_1_140421 | ERR13120660 |
| A2SW061221 | A_2_061221 | ERR17468632 |
| 1FW021221 | D_1_021221 | ERR17465016 |

Freshwater/FW is the freshwater site D in the final manuscript. A/B/C retain
their site letters; 1/2/3 retain their filter number; six-digit dates retain
their DDMMYY representation. Numeric ENA sample aliases were not used.
`bin/map_ena_titles.py` implements these explicit rules and checks a bijection
against all 60 canonical filters, all 20 biological samples and all dates.

To regenerate the crosswalk from the saved report:

```bash
python bin/map_ena_titles.py --filters metadata/filters.tsv \
  --ena-report metadata/ena_PRJEB75359_2026-11-12.tsv --output filters.remapped.tsv
```

Runtime resolution rechecks study membership and FASTQ checksums. To use the
snapshot explicitly, supply `--ena_report metadata/ena_PRJEB75359_2026-11-12.tsv`
to Nextflow or `--ena-report` to the Python resolver. The ENA report is metadata,
not a download of the sequencing reads.
