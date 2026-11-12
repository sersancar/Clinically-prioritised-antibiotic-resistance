# Validation scope — 12 November 2026

The rebuilt workflow is tested against methods and synthetic truth, without
asserting historical result totals or retaining the former scripts.

Passed locally:

* Nextflow 25.04.8 full wiring run: nineteen task instances completed with stubs.
* Nine Python method tests, including the study/run crosswalk, atomic MD5
  rejection, custom Kraken scoring, pair-name handling, rank-aware taxonomy,
  pooling/denominators, missing-filter rejection, BAM overlap coverage and
  unverifiable SNP-target exclusion.
* A real AMR++ SNP verification test containing five mutant and five wild-type
  reads: confirmed count five from baseline ten. The pinned upstream scorer
  ran using the socket-free serial scheduler.
* A real synthetic paired-read tool chain: FastQC, fastp, MultiQC, BWA-MEM,
  SAMtools duplicate removal/counts/coverage, Bowtie2, Kraken2 with a newly built
  tiny database, taxonomy linkage, Bracken distribution construction and
  species estimation. Two qualifying targets and 122 unique query/reference
  associations per catalogue were recovered, with two species estimated.
* MEGARes 3.00 automatic retrieval and resource preparation with local CARPDM/
  SNP inputs, frozen taxonomy and the NC_024711.1 reference retrieved from NCBI.
* Complete synthetic 60-filter aggregation and R reporting with default 9999
  requested permutations. All eleven required PDFs and thirty-three TSV table
  components/additional summaries were generated. Tests checked all six saved
  permutation matrices remain within date, the expected model/contrast sets and
  the three-test environmental BH family.
* Both declared Conda specifications solve on linux-64. Explicit resolved
  package manifests are included under envs/. R packages and the read-processing
  tools were exercised locally; Kraken2 and Bracken were compiled from their
  official releases for the tool smoke test.
* A second full report with no qualifying CARPDM targets completed, retaining
  blank completeness percentages and informative empty screening panels.

The complete multi-terabyte ENA production dataset has not been executed here.
The wiring run uses stubs and is not a raw-read scientific validation run. Real
alignment/SNP/classification tests use deliberately small synthetic data. A
successful Conda solve does not prove every package archive can be downloaded
on every host. Keep runtime versions/checksums and inspect the first production
run before interpreting study-specific outcomes.
