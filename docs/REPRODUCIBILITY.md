# Scientific and release validation

## Implemented rules

* 60 filters, 20 site-date biological samples, three filters each. Site C has
  two winter dates; A/B/D have six matched dates. No geographic coordinates.
* ENA run selection is explicit, study membership checked, exactly two mates
  required per run, MD5 checked before files become available. Multiple runs
  for one filter are concatenated in consistent run order.
* fastp Q20, at most 20% unqualified bases, at most five Ns, length >=50.
  Runtime `summary.after_filtering.total_reads` counts read ends across mates.
* Bowtie2 alignments: MAPQ >=30, `-F 2308` excludes unmapped, secondary and
  supplementary records; CARPDM breadth >=0.8 and mean depth >=5. Relaxed
  screening omits only depth. No CARPDM duplicate removal is added.
  The supplied `CARDGenes.sh` additionally establishes end-to-end mapping,
  no mixed/discordant alignments, `-k 1` and proper-pair filtering (`-f 2`),
  now applied by the CARPDM process.
* CARPDM SAMtools depth has no `-s`, includes full reference denominators and
  no depth cap. `-q` is base quality and `-Q` mapping quality for **depth**.
  Marker **mpileup** instead uses `-q 30 -Q 0`, as in the supplement.
* Nonqualifying filters supply zero numerator but remain in each biological
  denominator. RPM uses pooled read ends; RPKM additionally divides by kb.
* Reference FASTA checksums are mandatory for real screening runs.

## Numerical release gates after integrating original code

The final supplied supplement gives the following checks, to be evaluated
against recomputed data, never inserted as generated output:

| Quantity | Primary | Relaxed |
|---|---:|---:|
| Unique CARPDM references | 67 | 110 |
| Filter-reference calls | 820 | 1,303 |
| Biological sample-reference detections | 324 | 554 |

Primary associations: 333,616 mapped read ends; 166,808 distinct query-name/
reference associations; 1,558 genus-resolved, 81,190 above genus, 17,015 root,
67,045 with no retained assignment. Sites A/C have no primary detections.
crAssphage high-coverage qualification: 18 B filters and no A/C/D filters.

Check every filter total and reference metric before comparisons of pooled
tables. Compare numerical tables by exact identifiers and documented tolerances,
not visual similarities. Permutation p values require the original seed,
permutation matrix and package versions. Reconcile discrepancies rather than
changing results to satisfy these numbers.

## Statistics to preserve in the original report code

Hellinger gene profiles; uncorrected Bray-Curtis for `adonis2(Date + Site,
by='margin')`; permutations blocked within dates; bias-adjusted `betadisper`
with the same scheme. Cailliez correction only for PCoA coordinates. Sampling
subsets: ABD six dates and ABCD two common winter dates, RPM and RPKM. No
dispersion inference for the two-sample groups.

Friedman tests blocked by date; exact paired mean-difference tests over all
sign changes; BH within outcome/subset. Separate partial RDAs of standardised
water temperature, pH and salinity conditioned on Site + Date, date-restricted
permutations, BH over those three tests. No independent rainfall/humidity or
other date-only tests after conditioning on date.

Taxonomic affiliations describe sequence context, not physical host linkage.
Environmental associations do not establish a climate-change effect.

## Current test scope

Python tests exercise design validation, sample mapping, paired ENA file checks,
atomic checksum validation, coverage denominators, zero references, threshold
rules, complete-filter requirements and pooling denominators. Nextflow
`-stub-run` tests screening process wiring with synthetic outputs. Neither test
establishes numerical manuscript reproduction or validates real tool commands.
The conda specification is transcribed from Table S3 and remains unresolved.
