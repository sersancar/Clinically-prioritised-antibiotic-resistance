# Implemented methods

## Units and normalisation

The study manifest has 60 within-sample filters and 20 biological site–date
samples. A/B/D span six dates; C spans two winter dates. Runtime post-QC read-end
totals come from fastp. Published read totals are not used as denominators.
Biological RPM is 10^6 times summed eligible counts divided by the read-end
total of **all three filters**. Failed filters contribute zero numerator, not a
smaller denominator. RPKM divides reference RPM by length in kb before gene or
class aggregation. Neither measure estimates cells, copies per litre or biomass.

Broad read-end counts use deduplicated primary MAPQ >=30 alignments and exact
breadth >=0.8. BWA nonprimary/unmapped records are discarded before sorting to
avoid carrying unmapped reads through large temporary BAMs. SAMtools markdup -r
removes duplicates. SNP-dependent targets are verified using the upstream
AMR++ SNP engine, pinned to commit
8401ce3fff2860b257a9bcef23ac155e1e3e5222. Only eligible positive targets are sent
to the verifier. Targets requiring SNP confirmation but lacking verification
information receive zero confirmed count. Verification failures stop execution.

The SNP engine's read classification/count adjustment code is unchanged. A
serial per-gene scheduler replaces its nested process pool; Nextflow already
parallelises libraries. Source files and reference checksums are recorded.
Confirmed SNP counts can represent the engine's adjusted read-end estimate.
Taxonomic affiliations include alignments to retained references; they are not
restricted to individually mutation-confirmed query identifiers.

CARPDM uses Bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -k 1,
proper-pair filtering, MAPQ >=30 and primary alignments. Breadth uses full FASTA
length; mean depth is summed per-base coverage/full length. SAMtools depth is
used with no overlap suppression. The primary rule is breadth >=0.8 and mean
depth >=5; relaxed screening omits only depth. Rules are applied per filter.
Exact breadth is used, with no historical percentage rounding before screening.
These are operational screens, not calibrated detection limits.

crAssphage uses the separate very-sensitive mapping specification, MAPQ >=30,
primary alignments and samtools mpileup -aa -B -d 0 -q 30 -Q 0. Marker RPM pools
all retained alignments, including low-coverage filters. Coverage qualification
is reported separately. No unsupported ARG/marker dilution comparison is made.

## Taxonomy

Kraken2 classifies paired reads jointly. The custom score is assigned-taxid
k-mer counts divided by counts assigned to nonzero nonambiguous taxids.
Calculable scores below 0.2 are rejected; uncalculable scores are retained.
This is not the Kraken2 confidence parameter. Reports reconstructed from retained
classifications feed Bracken at species rank, read length 150, threshold ten.
Community estimates count read pairs; community biological tables label those
units separately from resistance read-end counts.

Only ARG-associated query identifiers are stored in memory. A per-filter SQLite
join deduplicates mates into unique query/reference associations. The classifier
output is scanned once. Terminal /1 and /2 mate suffixes are normalised; arbitrary
read-name substrings are not rewritten. Different references remain different
associations. NCBI ancestry, rather than whitespace or lineage position, resolves
ranks and domains. A taxid absent from the matching frozen taxonomy stops the
run. Root-only and missing/rejected classifications remain in completeness
denominators. Percentages with no eligible associations are blank.

S3 uses the declared 57-genus candidate list, without abundance exclusion.
Networks use pooled association RPM. Persistent nodes are positive on every
sampled date, then retained only above the site/node-type median. Bacteria and
Archaea are the eligible taxonomic nodes. Resistance-class totals retain all
original partners; isolated retained nodes remain in the display. Edges need
only be detected on at least one date. Network geometry and colours may differ
from the original figure; all selection data are exported.

These affiliations do not identify genomic hosts, pathogenic strains, viable
organisms, plasmid/chromosomal carriage or horizontal transfer.

## Statistical models

All multivariate/alpha tests use biological samples. Gene profiles are
Hellinger-transformed; Bray–Curtis distances feed marginal site tests in the
Date + Site PERMANOVA model. Permutations are restricted within date. Requested
permutations default to 9999; complete enumeration is used when appropriate.
Permutation matrices are saved. Cailliez correction applies only to displayed
PCoA coordinates, not the tested distances. Descriptive normal data ellipses
are restricted to nondegenerate groups with at least four samples.

Dispersion uses bias-adjusted betadisper under the same blocked permutation
scheme. Two-sample-per-site winter subsets have no dispersion hypothesis test;
small groups in the all-sample result are explicitly flagged. Sampling and
normalisation sensitivities compare ABD over six dates and all four sites over
two common dates, using RPM and RPKM.

Alpha outcomes are summed abundance, richness, exponential Shannon and inverse
Simpson. Friedman tests block by date. Exact paired contrasts enumerate all
sign changes and compare absolute mean differences; BH applies within outcome,
normalisation and sampling subset. Sequencing-depth Spearman tests are descriptive.

Separate partial RDAs test standardised water temperature, pH and salinity after
conditioning on site and date, using date-restricted permutations. BH correction
covers those three tests. Shared date-level weather/tide predictors receive no
independent hypothesis tests after conditioning on date. Non-estimable predictors
are reported rather than forced into a model. Environmental associations do not
establish climate-change causation.

## Sources

* [AMR++](https://github.com/Microbial-Ecology-Group/AMRplusplus), MEGARes 3.00
  nucleotide reference from commit fded6074d493fcdaf1ffd0c192b6b794249a464d.
* [AMR++ SNP verification](https://github.com/Microbial-Ecology-Group/AmrPlusPlus_SNP).
* [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/) for study reads.
* [NCBI nuccore](https://www.ncbi.nlm.nih.gov/nuccore/NC_024711.1) for the marker.

CARPDM and Kraken2/frozen taxdump inputs are supplied biological resources,
recorded with runtime provenance. They are not recovered from historical scripts.
