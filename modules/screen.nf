process BUILD_INDEX {
    tag "${label}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(label), path(fasta)
    output:
    tuple val(label), path(fasta), path('index')
    script:
    """
    mkdir index
    bowtie2-build --threads ${task.cpus} '${fasta}' index/reference
    """
    stub:
    """
    mkdir index
    touch index/reference.1.bt2
    """
}

process MAP_SCREEN {
    tag "${filter_id}:${label}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/screening", mode: 'copy', pattern: '*.metrics.tsv'
    input:
    tuple val(filter_id), path(reads), val(label), path(fasta), path(index)
    output:
    tuple val(label), path("${filter_id}.${label}.metrics.tsv"), emit: metrics
    tuple val(filter_id), val(label), path("${filter_id}.${label}.bam"), emit: bam
    path '*.bowtie2.log', emit: logs
    script:
    def coverage = label == 'crass' ?
        "samtools mpileup -aa -B -d 0 -q 30 -Q 0 '${filter_id}.${label}.bam' > depth.tsv" :
        "samtools depth -aa -d 0 -q 0 -Q 30 '${filter_id}.${label}.bam' > depth.tsv"
    def pileupFlag = label == 'crass' ? '--pileup' : ''
    def pairedOptions = label == 'carpdm' ? '--end-to-end --no-mixed --no-discordant -k 1' : ''
    def pairFilter = label == 'carpdm' ? '-f 2' : ''
    """
    bowtie2 --${params.bowtie2_preset} ${pairedOptions} -p ${task.cpus} -x '${index}/reference' \
      -1 '${reads[0]}' -2 '${reads[1]}' 2> ${filter_id}.${label}.bowtie2.log | \
      samtools view -b ${pairFilter} -q 30 -F 2308 - | \
      samtools sort -@ ${task.cpus} -o ${filter_id}.${label}.bam -
    ${coverage}
    samtools view ${filter_id}.${label}.bam > alignments.sam
    coverage_metrics.py --fasta '${fasta}' --depth depth.tsv --sam alignments.sam \
      --filter-id '${filter_id}' ${pileupFlag} --output ${filter_id}.${label}.metrics.tsv
    rm alignments.sam depth.tsv
    """
    stub:
    """
    printf 'filter_id\treference\tlength\tcovered_bases\tbreadth\tmean_depth\tmapped_read_ends\tprimary\trelaxed\n${filter_id}\ttest\t100\t0\t0\t0\t0\tFalse\tFalse\n' > ${filter_id}.${label}.metrics.tsv
    touch ${filter_id}.${label}.bam ${filter_id}.${label}.bowtie2.log
    """
}

process POOL_SCREEN {
    tag "${label}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/screening", mode: 'copy'
    input:
    tuple val(label), path(metrics)
    path json
    path filters
    output:
    tuple val(label), path("${label}.biological.tsv")
    script:
    def m = metrics.collect { "'${it}'" }.join(' ')
    def j = json.collect { "'${it}'" }.join(' ')
    """
    pool_screening.py --filters '${filters}' --metrics ${m} --fastp ${j} --output ${label}.biological.tsv
    """
    stub:
    """
    echo 'STUB ONLY: no biological results' > ${label}.biological.tsv
    """
}
