process VALIDATE_AMR_MANIFEST {
    input:
    path filters, stageAs: 'filter_metadata.tsv'
    path manifest, stageAs: 'bam_manifest.tsv'
    output:
    path 'validated_bams.tsv'
    script:
    """
    validate_amr_manifest.py --filters '${filters}' --manifest '${manifest}' --output validated_bams.tsv
    """
    stub:
    """
    printf 'filter_id\tbam\nA_1_290321\t${projectDir}/tests/reference.fa\n' > validated_bams.tsv
    """
}

process REBUILD_AMR_COUNTS {
    tag "${filter_id}:${params.amr_count_policy}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/amr_counts", mode: 'copy', pattern: '*.amr/*'
    input:
    tuple val(filter_id), path(bam)
    output:
    tuple val(filter_id), path("${filter_id}.amr"), emit: products
    path '*.amr/*.tsv', emit: tables
    path '*.amr/provenance.json', emit: provenance
    path '*.amr/samtools_version.txt', emit: versions
    script:
    def exclude = params.amr_count_policy == 'legacy2308' ? 2308 : 3332
    def depthExclude = params.amr_count_policy == 'legacy2308' ? 2048 : 3332
    """
    mkdir ${filter_id}.amr
    samtools depth -@ ${task.cpus} -aa -Q 30 -G ${depthExclude} '${bam}' > depth.tsv
    samtools view -F ${exclude} -q 30 '${bam}' | \
      count_amr_bam.py --depth depth.tsv --sam - --policy '${params.amr_count_policy}' --output-dir ${filter_id}.amr
    samtools --version > ${filter_id}.amr/samtools_version.txt
    ln -s "\$(realpath '${bam}')" ${filter_id}.amr/deduplicated.bam
    rm depth.tsv
    """
    stub:
    """
    mkdir ${filter_id}.amr
    printf 'gene_accession\tcount\n' > ${filter_id}.amr/reference_counts.tsv
    echo 'STUB ONLY' > ${filter_id}.amr/reference_metrics.tsv
    echo '{"stub":true}' > ${filter_id}.amr/provenance.json
    echo 'STUB ONLY' > ${filter_id}.amr/samtools_version.txt
    """
}
