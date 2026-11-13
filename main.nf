nextflow.enable.dsl=2
include { RESOLVE_ENA; DOWNLOAD; MERGE_READS; FASTQC; FASTP; MULTIQC } from './modules/reads'
include { PREPARE_KRAKEN; PREPARE_RESOURCES; BUILD_BWA; BUILD_BOWTIE; MAP_BROAD; VERIFY_SNPS; MAP_SCREEN; KRAKEN; LINK_READS; BRACKEN; ASSEMBLE_TABLES; REPORT } from './modules/analysis'

workflow {
    if (!(params.mode in ['full','download'])) error 'mode must be full or download'
    def filters=file(params.filters,checkIfExists:true)
    def ready=Channel.value(true)
    def resources=Channel.empty()
    def database=null
    def taxonomy=null
    def krakenProvenance=Channel.empty()
    if(params.mode=='full') {
        if(!workflow.stubRun && !params.carpdm_fasta) error 'Supply --carpdm_fasta or use data/carpdm_genes.fna.gz'
        def carp=file(workflow.stubRun ? "${projectDir}/tests/reference.fa" : params.carpdm_fasta,checkIfExists:true)
        if(params.kraken_db) {
            database=Channel.value(file(params.kraken_db,checkIfExists:true))
            def localDb=file(params.kraken_db,checkIfExists:true)
            def taxPath=params.taxonomy_dir ?: (file("${localDb}/nodes.dmp").exists() ? localDb : "${localDb}/taxonomy")
            taxonomy=Channel.value(file(workflow.stubRun ? "${projectDir}/tests" : taxPath,checkIfExists:true))
            if(!workflow.stubRun) {
                ['hash.k2d','opts.k2d','taxo.k2d','database150mers.kmer_distrib'].each { name ->
                    if(!file("${localDb}/${name}").exists()) error "Kraken2/Bracken database missing ${name}"
                }
            }
        } else {
            if(params.taxonomy_dir) error '--taxonomy_dir requires --kraken_db; automatic download includes matching taxonomy'
            PREPARE_KRAKEN(params.kraken_url,params.kraken_checksums_url,params.kraken_manifest_sha256)
            krakenProvenance=PREPARE_KRAKEN.out.provenance
            database=PREPARE_KRAKEN.out.database
            taxonomy=PREPARE_KRAKEN.out.taxonomy
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
            VERIFY_SNPS.out.status,LINK_READS.out.links,BRACKEN.out.species,BRACKEN.out.versions,PREPARE_RESOURCES.out,krakenProvenance).collect()
        ASSEMBLE_TABLES(artifacts,filters,file("${projectDir}/data",checkIfExists:true))
        REPORT(ASSEMBLE_TABLES.out,file("${projectDir}/analysis/report.R"))
    }
}
