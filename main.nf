nextflow.enable.dsl=2
include { RESOLVE_ENA; DOWNLOAD; MERGE_READS; FASTQC; FASTP; MULTIQC } from './modules/reads'
include { PREPARE_RESOURCES; BUILD_BWA; BUILD_BOWTIE; MAP_BROAD; VERIFY_SNPS; MAP_SCREEN; KRAKEN; LINK_READS; BRACKEN; ASSEMBLE_TABLES; REPORT } from './modules/analysis'

workflow {
    if (!(params.mode in ['full','download'])) error 'mode must be full or download'
    def filters=file(params.filters,checkIfExists:true)
    def ready=Channel.value(true)
    def resources=Channel.empty()
    def database=null
    def taxonomy=null
    if(params.mode=='full') {
        if(!workflow.stubRun && (!params.carpdm_fasta || !params.kraken_db)) {
            error 'Supply --carpdm_fasta and --kraken_db; no historical scripts/configurations are required'
        }
        def carp=file(params.carpdm_fasta ?: "${projectDir}/tests/reference.fa",checkIfExists:true)
        database=file(params.kraken_db ?: "${projectDir}/tests",checkIfExists:true)
        def taxPath=params.taxonomy_dir ?: (file("${database}/nodes.dmp").exists() ? database : "${database}/taxonomy")
        taxonomy=file(workflow.stubRun ? "${projectDir}/tests" : taxPath,checkIfExists:true)
        if(!workflow.stubRun) {
            ['hash.k2d','opts.k2d','taxo.k2d','database150mers.kmer_distrib'].each { name ->
                if(!file("${database}/${name}").exists()) error "Kraken2/Bracken database missing ${name}"
            }
        }
        def absent=file("${projectDir}/assets/NO_INPUT")
        PREPARE_RESOURCES(carp,params.megares_fasta?file(params.megares_fasta,checkIfExists:true):absent,
            params.crass_fasta?file(params.crass_fasta,checkIfExists:true):absent,
            params.snp_bundle?file(params.snp_bundle,checkIfExists:true):absent,taxonomy)
        resources=PREPARE_RESOURCES.out
        ready=resources.map{ true }
        BUILD_BWA(resources)
        BUILD_BOWTIE(PREPARE_RESOURCES.out.flatMap { r -> [tuple('carpdm',r.resolve('carpdm.fa')),tuple('crass',r.resolve('crass.fa'))] })
    }
    RESOLVE_ENA(filters,ready)
    downloads=RESOLVE_ENA.out.downloads.splitCsv(header:true,sep:'\t').map { r ->
        tuple(r.filter_id,r.run_accession,r.url1,r.url2,r.md5_1,r.md5_2)
    }
    DOWNLOAD(downloads)
    MERGE_READS(DOWNLOAD.out.groupTuple())
    if(params.mode=='download') {
        MERGE_READS.out.view { id,reads -> "Verified paired FASTQs for ${id}: ${reads}" }
    } else {
        FASTQC(MERGE_READS.out)
        FASTP(MERGE_READS.out)
        MULTIQC(FASTQC.out.mix(FASTP.out.json,FASTP.out.html).collect())
        MAP_BROAD(FASTP.out.reads.combine(BUILD_BWA.out),resources)
        VERIFY_SNPS(MAP_BROAD.out.samples,resources)
        MAP_SCREEN(FASTP.out.reads.combine(BUILD_BOWTIE.out))
        KRAKEN(FASTP.out.reads,database)
        carp=MAP_SCREEN.out.samples.filter { id,label,bam,metrics -> label=='carpdm' }.map { id,label,bam,metrics -> tuple(id,bam,metrics) }
        LINK_READS(VERIFY_SNPS.out.samples.join(carp).join(KRAKEN.out.classifications),taxonomy)
        BRACKEN(LINK_READS.out.community,database)
        artifacts=FASTP.out.json.mix(FASTP.out.versions,MAP_BROAD.out.versions,KRAKEN.out.versions,
            MAP_SCREEN.out.samples.map { id,label,bam,metrics -> metrics },VERIFY_SNPS.out.counts,
            VERIFY_SNPS.out.status,LINK_READS.out.links,BRACKEN.out.species,BRACKEN.out.versions,PREPARE_RESOURCES.out).collect()
        ASSEMBLE_TABLES(artifacts,filters,file("${projectDir}/metadata"))
        REPORT(ASSEMBLE_TABLES.out,file("${projectDir}/analysis/report.R"))
    }
}
