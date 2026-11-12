// These adapters must be implemented from the author's original scripts.
// No generic replacement is claimed to reproduce their results.
process AMR_ORIGINAL {
    tag "${filter_id}"
    input:
    tuple val(filter_id), path(reads)
    path adapter
    output:
    tuple val(filter_id), path("${filter_id}.amr")
    script:
    """
    mkdir ${filter_id}.amr
    bash '${adapter}' '${filter_id}' '${reads[0]}' '${reads[1]}' '${filter_id}.amr'
    test -s ${filter_id}.amr/reference_counts.tsv
    test -s ${filter_id}.amr/deduplicated.bam
    test -s ${filter_id}.amr/provenance.json
    """
    stub:
    """
    mkdir ${filter_id}.amr
    echo 'STUB ONLY' > ${filter_id}.amr/reference_counts.tsv
    echo 'STUB ONLY' > ${filter_id}.amr/deduplicated.bam
    echo '{"stub":true}' > ${filter_id}.amr/provenance.json
    """
}

process KRAKEN {
    tag "${filter_id}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(filter_id), path(reads)
    path db
    output:
    tuple val(filter_id), path("${filter_id}.kraken.tsv")
    script:
    """
    kraken2 --db '${db}' --paired --use-names --threads ${task.cpus} --gzip-compressed \
      --output ${filter_id}.kraken.tsv '${reads[0]}' '${reads[1]}'
    """
    stub:
    """
    echo 'STUB ONLY' > ${filter_id}.kraken.tsv
    """
}

process TAXONOMY_ORIGINAL {
    tag "${filter_id}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(filter_id), path(classifications)
    path db
    path adapter
    output:
    tuple val(filter_id), path("${filter_id}.taxonomy")
    script:
    """
    mkdir ${filter_id}.taxonomy
    bash '${adapter}' '${filter_id}' '${classifications}' '${db}' '${filter_id}.taxonomy'
    test -f ${filter_id}.taxonomy/retained_assignments.tsv
    test -s ${filter_id}.taxonomy/bracken_species.tsv
    test -s ${filter_id}.taxonomy/provenance.json
    """
    stub:
    """
    mkdir ${filter_id}.taxonomy
    echo 'STUB ONLY' > ${filter_id}.taxonomy/retained_assignments.tsv
    echo 'STUB ONLY' > ${filter_id}.taxonomy/bracken_species.tsv
    echo '{"stub":true}' > ${filter_id}.taxonomy/provenance.json
    """
}

process REPORT_ORIGINAL {
    input:
    path adapter
    path filters
    path amr, stageAs: 'amr/*'
    path taxonomy, stageAs: 'taxonomy/*'
    path bams, stageAs: 'screen_bams/*'
    path metrics, stageAs: 'screen_metrics/*'
    path biological, stageAs: 'biological/*'
    path qc, stageAs: 'qc/*'
    path provenance
    output:
    path 'manuscript'
    publishDir "${params.outdir}", mode: 'copy'
    script:
    """
    mkdir manuscript
    bash '${adapter}' '${filters}' amr taxonomy screen_bams screen_metrics biological qc manuscript
    validate_outputs.py manuscript '${projectDir}/assets/expected_outputs.json'
    cp '${provenance}' manuscript/reference_provenance.json
    """
    stub:
    """
    mkdir manuscript
    echo 'Synthetic wiring check only; no manuscript outputs generated.' > manuscript/STUB_ONLY.txt
    """
}
