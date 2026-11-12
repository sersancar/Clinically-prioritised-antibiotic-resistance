# Development validation — 12 November 2026

Passed locally:

* Eighteen Python method/integration tests (`python -m unittest discover -s tests -v`).
* Python compilation of `bin/` and `tests/`.
* Full Nextflow synthetic wiring run under **25.04.8**, the version recorded
  in Table S3: all 17 process tasks completed.
* Full synthetic wiring run under Nextflow 25.10.0, plus a separate screening
  synthetic run after adding the reference-verification dependency.
* Live ENA metadata retrieval for PRJEB75359 and bijective matching of all
  60 run sample titles to all 60 canonical filters.
* Offline resolution of the saved ENA report: all runs paired, mate order and
  checksum format validated. File checksums are checked only on actual download.

Not executed or validated:

* Production sequencing downloads, alignment, taxonomic classification, real
  SAMtools depth/mpileup compatibility or original AMR++ execution.
* Conda solving, container builds, Slurm execution or performance estimates.
* Original R analyses, figures, Cytoscape layouts or final Excel packaging.
* Numerical agreement with published screening, taxonomy or statistical results.
* Remote GitHub push (no authenticated GitHub session).

Synthetic wiring checks deliberately replace scientific processes with fixtures
and generate no manuscript figures or tables. They validate dataflow and file
staging, not scientific correctness. See `OUTPUT_STATUS.md` for remaining work.

## Uploaded-script integration checks

Passed: scoring equivalence against the original embedded Python block; report
clade counts and names on a synthetic taxonomy; source archive SHA256 checks;
60-filter R-input conversion using synthetic runtime totals. All uploaded Bash
scripts and the new taxonomy adapter passed bash -n. The FROM_TABLES workflow
completed a three-task synthetic wiring run under Nextflow 25.04.8.

R and the production bioinformatics tools are unavailable in this workspace.
The portable R changes have not been executed or parsed with R locally. A CI
R syntax job was added but has not been run on GitHub. Dependency specifications
remain unresolved; no claim of numerical manuscript reproduction is made.

## Additional AMR evidence checks

Four more tests pass: launch argument preservation, duplicate-flag policy
differences, original AWK breadth-rounding equivalence and complete unique BAM
manifest validation. No real AMR++ or BAM counting was executed. The newly
received launch script confirms requested flags but not SNP effects on the
rebuilt matrix. See AMR_EXECUTION.md.

Nextflow 25.04.8 also completed the two-task FROM_AMR_BAMS synthetic wiring
check and the main full-mode synthetic workflow after these additions. These
checks execute fixture stubs, not production alignments/count reconstruction.

Resumed-run report filename collisions found during testing were addressed with
timestamped execution/trace/timeline/DAG filenames, preserving each run report.
