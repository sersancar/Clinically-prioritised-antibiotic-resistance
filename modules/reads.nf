process RESOLVE_ENA {
    label 'small'
    publishDir "${params.outdir}/provenance/ena", mode: 'copy'
    input:
    path filters
    val ready
    output:
    path 'ena/downloads.tsv', emit: downloads
    path 'ena/*', emit: provenance
    script:
    """
    ena_manifest.py --filters '${filters}' --study '${params.study}' --outdir ena
    """
    stub:
    """
    mkdir ena
    printf 'filter_id\trun_accession\turl1\turl2\tmd5_1\tmd5_2\nA_1_290321\tERR000001\thttps://ftp.sra.ebi.ac.uk/test_1.fastq.gz\thttps://ftp.sra.ebi.ac.uk/test_2.fastq.gz\t00000000000000000000000000000000\t00000000000000000000000000000000\n' > ena/downloads.tsv
    """
}
process DOWNLOAD {
    tag "${id}:${run}"
    label 'small'
    maxForks 4
    errorStrategy 'retry'
    maxRetries 2
    input:
    tuple val(id),val(run),val(url1),val(url2),val(md51),val(md52)
    output:
    tuple val(id),path("${run}_1.fastq.gz"),path("${run}_2.fastq.gz")
    script:
    """
    download_fastq.py '${url1}' '${md51}' '${run}_1.fastq.gz'
    download_fastq.py '${url2}' '${md52}' '${run}_2.fastq.gz'
    """
    stub:
    """
    printf '@test\nACGT\n+\nIIII\n' | gzip > ${run}_1.fastq.gz
    cp ${run}_1.fastq.gz ${run}_2.fastq.gz
    """
}
process MERGE_READS {
    tag "${id}"
    label 'small'
    input:
    tuple val(id),path(mates1),path(mates2)
    output:
    tuple val(id),path("${id}_R{1,2}.fastq.gz")
    script:
    def r1=mates1.sort{it.name}.collect{"'${it}'"}.join(' ')
    def r2=mates2.sort{it.name}.collect{"'${it}'"}.join(' ')
    def command1=mates1.size()==1?"ln -s \"\$(realpath '${mates1[0]}')\" ${id}_R1.fastq.gz":"cat ${r1} > ${id}_R1.fastq.gz"
    def command2=mates2.size()==1?"ln -s \"\$(realpath '${mates2[0]}')\" ${id}_R2.fastq.gz":"cat ${r2} > ${id}_R2.fastq.gz"
    """
    ${command1}
    ${command2}
    """
}
process FASTQC {
    tag "${id}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/qc/fastqc",mode:'copy'
    input:
    tuple val(id),path(reads)
    output:
    path '*_fastqc*'
    script:
    """
    fastqc -t ${task.cpus} '${reads[0]}' '${reads[1]}'
    """
    stub:
    """
    touch ${id}_fastqc.html ${id}_fastqc.zip
    """
}
process FASTP {
    tag "${id}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/qc/fastp",mode:'copy',pattern:'*.fastp.*'
    input:
    tuple val(id),path(reads)
    output:
    tuple val(id),path("${id}.clean_R{1,2}.fastq.gz"),emit:reads
    path '*.fastp.json',emit:json
    path '*.fastp.html',emit:html
    path '*.versions.json',emit:versions
    script:
    """
    fastp -i '${reads[0]}' -I '${reads[1]}' -o ${id}.clean_R1.fastq.gz -O ${id}.clean_R2.fastq.gz \
      --qualified_quality_phred 20 --unqualified_percent_limit 20 --n_base_limit 5 --length_required 50 \
      --thread ${task.cpus} --json ${id}.fastp.json --html ${id}.fastp.html
    software_versions.py --output ${id}.fastp.versions.json fastp
    """
    stub:
    """
    cp '${reads[0]}' ${id}.clean_R1.fastq.gz
    cp '${reads[1]}' ${id}.clean_R2.fastq.gz
    echo '{"summary":{"after_filtering":{"total_reads":2}}}' > ${id}.fastp.json
    touch ${id}.fastp.html
    echo '{"fastp":"STUB"}' > ${id}.fastp.versions.json
    """
}
process MULTIQC {
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/qc/multiqc",mode:'copy'
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
