process PREPARE_KRAKEN {
    label 'small'
    time '168h'
    input:
    val archive_url
    val checksums_url
    val manifest_sha256
    output:
    path 'kraken_database', emit: database
    path 'kraken_taxonomy', emit: taxonomy
    path 'kraken_database.json', emit: provenance
    script:
    """
    prepare_kraken.py --url '${archive_url}' --checksums-url '${checksums_url}' \
      --manifest-sha256 '${manifest_sha256}' --database kraken_database \
      --taxonomy kraken_taxonomy --provenance kraken_database.json
    """
    stub:
    """
    mkdir kraken_database kraken_taxonomy
    touch kraken_database/hash.k2d kraken_database/opts.k2d kraken_database/taxo.k2d
    touch kraken_database/database150mers.kmer_distrib kraken_taxonomy/nodes.dmp kraken_taxonomy/names.dmp
    echo '{"stub":true}' > kraken_database.json
    """
}

process PREPARE_RESOURCES {
    label 'small'
    input:
    path carpdm,stageAs:'input_carpdm.fa'
    path megares,stageAs:'input_megares'
    path crass,stageAs:'input_crass'
    path snp,stageAs:'input_snp'
    path taxonomy
    output:
    path 'resources'
    script:
    // Optional inputs are staged with fixed names; pass sentinel status explicitly.
    def megArg=params.megares_fasta?megares:"${projectDir}/assets/NO_INPUT"
    def crArg=params.crass_fasta?crass:"${projectDir}/assets/NO_INPUT"
    def snpArg=params.snp_bundle?snp:"${projectDir}/assets/NO_INPUT"
    """
    resources.py --carpdm '${carpdm}' --megares '${megArg}' --crass '${crArg}' --snp-bundle '${snpArg}' \
      --taxonomy '${taxonomy}' --out resources
    mv resources/provenance.json resources/resource_provenance.json
    """
    stub:
    """
    mkdir -p resources/snp
    cp '${projectDir}/tests/reference.fa' resources/carpdm.fa
    cp '${projectDir}/tests/reference.fa' resources/crass.fa
    cp '${projectDir}/tests/reference.fa' resources/megares.fa
    echo '{"stub":true}' > resources/resource_provenance.json
    """
}
process BUILD_BWA {
    conda "${projectDir}/envs/reads.yml"
    input:
    path resources
    output:
    path 'bwa_index'
    script:
    """
    mkdir bwa_index
    cp '${resources}/megares.fa' bwa_index/megares.fa
    bwa index bwa_index/megares.fa
    """
    stub:
    """
    mkdir bwa_index
    touch bwa_index/megares.fa
    """
}
process BUILD_BOWTIE {
    tag "${catalogue}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(catalogue),path(fasta)
    output:
    tuple val(catalogue),path(fasta),path('index')
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
process MAP_BROAD {
    tag "${id}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(id),path(reads),path(index)
    path resources
    output:
    tuple val(id),path("${id}.broad.bam"),path("${id}.broad.metrics.tsv"),emit:samples
    path '*.versions.json',emit:versions
    script:
    """
    bwa mem -t ${task.cpus} -K 100000000 '${index}/megares.fa' '${reads[0]}' '${reads[1]}' | \
      samtools view -u -F 2308 - | \
      samtools sort -n -@ ${task.cpus} -o names.bam -
    samtools fixmate -m names.bam fixed.bam
    rm names.bam
    samtools sort -@ ${task.cpus} -o coordinate.bam fixed.bam
    rm fixed.bam
    samtools markdup -r coordinate.bam deduplicated.bam
    rm coordinate.bam
    samtools view -b -F 2308 -q 30 deduplicated.bam > ${id}.broad.bam
    rm deduplicated.bam
    samtools index ${id}.broad.bam
    bam_metrics.py --bam ${id}.broad.bam --fasta '${resources}/megares.fa' --filter-id '${id}' --catalogue broad --output ${id}.broad.metrics.tsv
    software_versions.py --output ${id}.broad.versions.json bwa samtools
    """
    stub:
    """
    touch ${id}.broad.bam ${id}.broad.metrics.tsv
    echo '{"samtools":"STUB"}' > ${id}.broad.versions.json
    """
}
process VERIFY_SNPS {
    tag "${id}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(id),path(bam),path(metrics)
    path resources
    output:
    tuple val(id),path("${id}.broad.bam"),path("${id}.broad.counts.tsv"),emit:samples
    path '*.broad.counts.tsv',emit:counts
    path '*.snp_status.json',emit:status
    script:
    """
    verify_snps.py --bam '${bam}' --metrics '${metrics}' --bundle '${resources}/snp' --fasta '${resources}/megares.fa' \
      --filter-id '${id}' --threads ${task.cpus} --out snp_check
    cp snp_check/counts.tsv ${id}.broad.counts.tsv
    cp snp_check/status.json ${id}.snp_status.json
    """
    stub:
    """
    touch ${id}.broad.counts.tsv
    echo '{"stub":true}' > ${id}.snp_status.json
    """
}
process MAP_SCREEN {
    tag "${id}:${catalogue}"
    conda "${projectDir}/envs/reads.yml"
    publishDir "${params.outdir}/screening",mode:'copy',pattern:'*.metrics.tsv'
    input:
    tuple val(id),path(reads),val(catalogue),path(fasta),path(index)
    output:
    tuple val(id),val(catalogue),path("${id}.${catalogue}.bam"),path("${id}.${catalogue}.metrics.tsv"),emit:samples
    script:
    def pairOptions=catalogue=='carpdm'?'--end-to-end --no-mixed --no-discordant -k 1':''
    def pairFilter=catalogue=='carpdm'?'-f 2':''
    """
    bowtie2 --very-sensitive ${pairOptions} -p ${task.cpus} -x '${index}/reference' -1 '${reads[0]}' -2 '${reads[1]}' | \
      samtools view -b ${pairFilter} -q 30 -F 2308 - | samtools sort -@ ${task.cpus} -o ${id}.${catalogue}.bam -
    samtools index ${id}.${catalogue}.bam
    bam_metrics.py --bam ${id}.${catalogue}.bam --fasta '${fasta}' --filter-id '${id}' --catalogue '${catalogue}' --output ${id}.${catalogue}.metrics.tsv
    """
    stub:
    """
    touch ${id}.${catalogue}.bam ${id}.${catalogue}.metrics.tsv
    """
}
process KRAKEN {
    tag "${id}"
    label 'taxonomy'
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(id),path(reads)
    path database
    output:
    tuple val(id),path("${id}.kraken.tsv"),emit:classifications
    path '*.versions.json',emit:versions
    script:
    def mappingValue=params.kraken_memory_mapping.toString().toLowerCase()
    if (!(mappingValue in ['true','false'])) error '--kraken_memory_mapping must be true or false'
    def mappingArg=mappingValue=='true' ? '--memory-mapping' : ''
    """
    kraken2 --db '${database}' ${mappingArg} --paired --gzip-compressed --threads ${task.cpus} \
      --output ${id}.kraken.tsv '${reads[0]}' '${reads[1]}'
    software_versions.py --output ${id}.taxonomy.versions.json kraken2
    """
    stub:
    """
    touch ${id}.kraken.tsv
    echo '{"kraken2":"STUB"}' > ${id}.taxonomy.versions.json
    """
}
process LINK_READS {
    tag "${id}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(id),path(broad_bam),path(counts),path(carp_bam),path(carp_metrics),path(kraken)
    path taxonomy
    output:
    path "${id}.links",emit:links
    tuple val(id),path("${id}.links/community.report"),emit:community
    script:
    """
    link_reads.py --filter-id '${id}' --broad-bam '${broad_bam}' --carpdm-bam '${carp_bam}' --broad-counts '${counts}' \
      --carpdm-metrics '${carp_metrics}' --kraken '${kraken}' --taxonomy '${taxonomy}' --out ${id}.links
    """
    stub:
    """
    mkdir ${id}.links
    touch ${id}.links/community.report
    """
}
process BRACKEN {
    tag "${id}"
    conda "${projectDir}/envs/reads.yml"
    input:
    tuple val(id),path(report)
    path database
    output:
    path '*.bracken.tsv',emit:species
    path '*.versions.json',emit:versions
    script:
    """
    if [ -s '${report}' ]; then
      bracken -d '${database}' -i '${report}' -o ${id}.bracken.tsv -r 150 -l S -t 10
    else
      printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n' > ${id}.bracken.tsv
    fi
    software_versions.py --output ${id}.bracken.versions.json bracken
    """
    stub:
    """
    touch ${id}.bracken.tsv
    echo '{"bracken":"STUB"}' > ${id}.bracken.versions.json
    """
}
process ASSEMBLE_TABLES {
    conda "${projectDir}/envs/reads.yml"
    input:
    path artifacts,stageAs:'inputs/*'
    path filters
    path metadata
    output:
    path 'report_data'
    script:
    """
    software_versions.py --output inputs/report.versions.json fastp fastqc multiqc bwa samtools bowtie2 kraken2 bracken
    aggregate_results.py --filters '${filters}' --inputs inputs --metadata '${metadata}' --out report_data
    """
    stub:
    """
    mkdir -p report_data/tables report_data/provenance
    """
}
process REPORT {
    conda "${projectDir}/envs/analysis.yml"
    publishDir "${params.outdir}",mode:'copy'
    input:
    path data
    path report_script
    output:
    path 'manuscript'
    script:
    """
    Rscript '${report_script}' '${data}' manuscript ${params.seed} ${params.permutations}
    cp -r '${data}/provenance' manuscript/
    validate_report.py manuscript
    """
    stub:
    """
    mkdir manuscript
    echo 'Synthetic wiring check only. No scientific results generated.' > manuscript/STUB_ONLY.txt
    """
}
