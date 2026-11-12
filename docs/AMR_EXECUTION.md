# AMR++ launcher and matrix-policy evidence

The additional uploads are archived unchanged as:

* `AMRplusplus_launch_ncahill.sh` (uploaded `AMR++(1).sh`);
* `CARDGenes_ncahill.sh` (uploaded `CARDGenes(1).sh`);
* `build_dedup_AMR_matrix_ncahill.sh`.

The source archive manifest now covers all fourteen files and their checksums.

## Confirmed original launcher

The launcher ran AMR++ with `--snp Y --deduped Y --pipeline resistome`, profile
`local_slurm`, paired `*_clean_{1,2}.fq.gz` reads, `-resume` and 36 allocated
threads. It selected Java from the AMR++ environment and preserved Nextflow work
files. It does not identify the AMR++ Git commit or archive the effective
configuration/resources. Those remain required for exact reproduction.

`bin/run_amrplusplus.py` provides an explicit external launch stage using those
scientific switches. It checks a full source commit and clean tracked source,
records source/additional configuration hashes and preserves the argument
vector as JSON. It avoids launching sixty hidden child Nextflow workflows.

```bash
python bin/run_amrplusplus.py \
  --source /absolute/path/AMRplusplus \
  --revision FULL_40_CHARACTER_GIT_COMMIT \
  --config /absolute/path/validated_amr.config \
  --reads-glob '/absolute/path/qc/*_clean_{1,2}.fq.gz' \
  --output /absolute/path/AMRpp_output \
  --profile local_slurm --threads 36 \
  --run-dir /absolute/path/AMRpp_execution --dry-run
```

Remove `--dry-run` only for an actual run on the intended host. Supply the
validated Nextflow/Java environment and a configuration with absolute resource
paths: the original launch directory was the AMR++ source checkout, whereas
this wrapper keeps its logs/work in the supplied execution directory. A matching
reference/SNP manifest is still required; configuration hashes alone do not
identify database contents. The pipeline's QC output naming must be exported
to the original `*_clean_1.fq.gz` / `*_clean_2.fq.gz` convention or the explicit
read glob must be adapted to the original AMR++ paired-read matcher.

## Two different count policies

| Policy | Original count exclusions (`samtools view -F`) | Original depth exclusions added with `-G` |
|---|---:|---:|
| `legacy2308` | 2308: unmapped, secondary, supplementary | 2048, in addition to SAMtools defaults |
| `exclude_duplicates3332` | 3332: the above plus duplicate-flagged records | 3332, in addition to defaults |

Both require MAPQ >=30. Both first print breadth as a percentage rounded to two
decimals, then apply >=80 to that rounded value. The portable counter preserves
this order. In a synthetic test, 16000/20001 covered bases is below 80% exactly
but rounds to 80.00%, so it passes the original rule. This boundary behaviour
must not be silently replaced by exact fraction screening when reproducing the
source matrix. The CARPDM fraction threshold is a separate rule.

No count policy is selected by default. The legacy depth command already
excludes duplicate-flagged records through SAMtools defaults while its read
counter can include them. The newer script excludes them explicitly from both.
The source BAMs' actual duplicate flags determine whether this changes results.

## Nextflow count rebuilding

Create a TSV with `filter_id,bam` columns for the sixty canonical manuscript
filters. Relative BAM paths resolve against the manifest's original directory;
all files must exist and no BAM may be assigned twice.

```bash
nextflow run main.nf -entry FROM_AMR_BAMS -profile conda,slurm \
  --amr_bam_manifest /absolute/path/deduplicated_bams.tsv \
  --amr_count_policy exclude_duplicates3332 \
  --outdir results/amr_rebuilt -resume
```

Use `legacy2308` instead to reproduce the earlier counting script. Each filter
emits reference counts, rounded breadth metrics, policy provenance and observed
SAMtools version. Counts are streamed from SAMtools rather than dumping all
alignments to a temporary SAM. Nextflow synchronises tasks instead of the source
script's last-array-task polling loop. Original BAMs remain unchanged.

The count-rebuild stage does not itself apply a SNP gate. The launcher confirms
that SNP screening was requested upstream, but does not establish how that
screen affected counts reconstructed from the retained BAM. Verify that link
before describing the final rebuilt matrix as SNP-confirmed.

## CARPDM and sample scope

The new CARDGenes version has the same mapping/threshold logic as the earlier
upload and activates an environment named `samtools_new`. That name does not
identify a version or prove whether its `coverage` or `depth` branch ran.
The different batch paths and array sizes (42 or 68) are not used to add samples
to this manuscript: all validation remains restricted to its 60-filter manifest.

Confirm which count-rebuild script produced the matrix underlying final Table
S8, and supply the effective AMR++ configuration/source commit/resource manifest.
These are the remaining AMR-specific decisions/evidence, rather than another
request for the launcher already supplied.
