process ORIGINAL_R_ANALYSIS {
    tag 'original analysis'
    conda "${projectDir}/envs/analysis.yml"
    publishDir "${params.outdir}/reanalysis", mode: 'copy'
    input:
    path inputs, stageAs: 'input'
    output:
    path 'revision'
    script:
    """
    export AREST_INPUT="\$PWD/input"
    export AREST_OUTPUT="\$PWD"
    export AREST_RUN="\$PWD/revision"
    export AREST_CARPDM_UNFILTERED='${params.carpdm_unfiltered ? 'TRUE' : 'FALSE'}'
    Rscript '${projectDir}/analysis/arestAnalysis_revised.R'
    validate_analysis_status.py revision/audit/module_status.tsv
    """
    stub:
    """
    mkdir -p revision/audit
    echo 'STUB ONLY: no R analysis performed' > revision/audit/STUB_ONLY.txt
    """
}

process ORIGINAL_R_FIGURES {
    conda "${projectDir}/envs/analysis.yml"
    publishDir "${params.outdir}/reanalysis", mode: 'copy'
    input:
    path revision
    output:
    path 'figures_original_numbering'
    script:
    """
    Rscript '${projectDir}/analysis/arest_revision_figures.R' '${revision}' figures_original_numbering
    """
    stub:
    """
    mkdir figures_original_numbering
    echo 'STUB ONLY: no figures generated' > figures_original_numbering/STUB_ONLY.txt
    """
}

process ORIGINAL_R_NETWORKS {
    conda "${projectDir}/envs/analysis.yml"
    publishDir "${params.outdir}/reanalysis", mode: 'copy'
    input:
    path revision
    path taxonomy, stageAs: 'broad_taxonomy'
    path domain_map, stageAs: 'domain_map.tsv'
    output:
    path 'networks_original_numbering'
    script:
    """
    cp -aL '${revision}' network_revision
    cp '${domain_map}' network_revision/taxon_domain_map.tsv
    export AREST_CYTOSCAPE=FALSE
    Rscript '${projectDir}/analysis/arest_networks.R' '${taxonomy}' network_revision networks_original_numbering
    """
    stub:
    """
    mkdir networks_original_numbering
    echo 'STUB ONLY: no networks generated' > networks_original_numbering/STUB_ONLY.txt
    """
}
