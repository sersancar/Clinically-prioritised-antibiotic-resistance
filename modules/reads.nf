process RESOLVE_ENA {
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/provenance/ena", mode: 'copy'
    input:
    path filters
    path report
    output:
    path 'ena/downloads.tsv', emit: downloads
    path 'ena/*', emit: provenance
    script:
    def reportArg = report.name == 'NO_ENA_REPORT' ? '' : "--ena-report '${report}'"
    """
    ena_manifest.py --filters '${filters}' --study '${params.study}' ${reportArg} --outdir ena
    """
    stub:
    """
    mkdir ena
    printf 'filter_id\trun_accession\turl1\turl2\tmd5_1\tmd5_2\nA_1_290321\tERR000001\thttps://ftp.sra.ebi.ac.uk/test_1.fastq.gz\thttps://ftp.sra.ebi.ac.uk/test_2.fastq.gz\t00000000000000000000000000000000\t00000000000000000000000000000000\n' > ena/downloads.tsv
    """
}

process DOWNLOAD {
    tag "${filter_id}:${run}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(filter_id), val(run), val(url1), val(url2), val(md5_1), val(md5_2)
    output:
    tuple val(filter_id), path("${run}_1.fastq.gz"), path("${run}_2.fastq.gz")
    script:
    """
    download_fastq.py '${url1}' '${md5_1}' '${run}_1.fastq.gz'
    download_fastq.py '${url2}' '${md5_2}' '${run}_2.fastq.gz'
    """
    stub:
    """
    printf '@read/1\nACGT\n+\nIIII\n' | gzip > ${run}_1.fastq.gz
    printf '@read/2\nACGT\n+\nIIII\n' | gzip > ${run}_2.fastq.gz
    """
}

process MERGE_READS {
    tag "${filter_id}"
    input:
    tuple val(filter_id), path(mates1), path(mates2)
    output:
    tuple val(filter_id), path("${filter_id}_R{1,2}.fastq.gz")
    script:
    // Lexicographic run order is identical for the two mates.
    def r1 = mates1.sort { it.name }.collect { "'${it}'" }.join(' ')
    def r2 = mates2.sort { it.name }.collect { "'${it}'" }.join(' ')
    """
    cat ${r1} > ${filter_id}_R1.fastq.gz
    cat ${r2} > ${filter_id}_R2.fastq.gz
    """
}

process FASTQC {
    tag "${filter_id}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    input:
    tuple val(filter_id), path(reads)
    output:
    path '*_fastqc*'
    script:
    """
    fastqc --threads ${task.cpus} '${reads[0]}' '${reads[1]}'
    """
    stub:
    """
    touch ${filter_id}_fastqc.html ${filter_id}_fastqc.zip
    """
}

process FASTP {
    tag "${filter_id}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.fastp.*'
    input:
    tuple val(filter_id), path(reads)
    output:
    tuple val(filter_id), path("${filter_id}.clean_R{1,2}.fastq.gz"), emit: reads
    path '*.fastp.json', emit: json
    path '*.fastp.html', emit: html
    script:
    """
    fastp --in1 '${reads[0]}' --in2 '${reads[1]}' \
      --out1 ${filter_id}.clean_R1.fastq.gz --out2 ${filter_id}.clean_R2.fastq.gz \
      --qualified_quality_phred 20 --unqualified_percent_limit 20 \
      --n_base_limit 5 --length_required 50 --thread ${task.cpus} \
      --json ${filter_id}.fastp.json --html ${filter_id}.fastp.html
    """
    stub:
    """
    cp '${reads[0]}' ${filter_id}.clean_R1.fastq.gz
    cp '${reads[1]}' ${filter_id}.clean_R2.fastq.gz
    echo '{"summary":{"after_filtering":{"total_reads":2}}}' > ${filter_id}.fastp.json
    touch ${filter_id}.fastp.html
    """
}

process MULTIQC {
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'
    input:
    path reports
    output:
    path 'multiqc*'
    script:
    """
    multiqc . --force
    """
    stub:
    """
    touch multiqc_report.html
    """
}
