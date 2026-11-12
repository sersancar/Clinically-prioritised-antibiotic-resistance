nextflow.enable.dsl=2
include { RESOLVE_ENA; DOWNLOAD; MERGE_READS; FASTQC; FASTP; MULTIQC } from './modules/reads'
include { BUILD_INDEX; MAP_SCREEN; POOL_SCREEN } from './modules/screen'
include { AMR_ORIGINAL; KRAKEN; TAXONOMY_ORIGINAL; REPORT_ORIGINAL } from './modules/adapters'
include { ORIGINAL_R_ANALYSIS; ORIGINAL_R_FIGURES; ORIGINAL_R_NETWORKS } from './modules/reanalysis'
include { VALIDATE_AMR_MANIFEST; REBUILD_AMR_COUNTS } from './modules/amr_counts'

workflow FROM_AMR_BAMS {
    if (!params.amr_bam_manifest) error 'Supply --amr_bam_manifest with filter_id,bam columns'
    if (!(params.amr_count_policy in ['legacy2308','exclude_duplicates3332'])) {
        error 'Select --amr_count_policy legacy2308 or exclude_duplicates3332; no scientific default is assumed'
    }
    VALIDATE_AMR_MANIFEST(file(params.filters,checkIfExists:true),file(params.amr_bam_manifest,checkIfExists:true))
    bams = VALIDATE_AMR_MANIFEST.out.splitCsv(header:true,sep:'\t').map { r ->
        tuple(r.filter_id,file(r.bam,checkIfExists:true))
    }
    REBUILD_AMR_COUNTS(bams)
}

workflow FROM_TABLES {
    if (!params.analysis_inputs) error 'FROM_TABLES requires --analysis_inputs; see docs/ORIGINAL_CODE_INTEGRATION.md'
    if ((params.broad_taxonomy && !params.domain_map) || (!params.broad_taxonomy && params.domain_map)) {
        error 'Supply --broad_taxonomy and --domain_map together for offline network reproduction'
    }
    ORIGINAL_R_ANALYSIS(file(params.analysis_inputs,checkIfExists:true))
    ORIGINAL_R_FIGURES(ORIGINAL_R_ANALYSIS.out)
    if (params.broad_taxonomy) {
        ORIGINAL_R_NETWORKS(ORIGINAL_R_ANALYSIS.out,
            file(params.broad_taxonomy,checkIfExists:true),file(params.domain_map,checkIfExists:true))
    }
}

process VERIFY_REFERENCES {
    input:
    path lock
    path carpdm, stageAs: 'carpdm.fa'
    path crass, stageAs: 'crass.fa'
    output:
    path 'reference_provenance.json'
    script:
    """
    verify_references.py --lock '${lock}' --carpdm '${carpdm}' --crass '${crass}' --output reference_provenance.json
    """
    stub:
    """
    echo '{"stub": true}' > reference_provenance.json
    """
}

workflow {
    if (!(params.mode in ['full','screening','download'])) error 'mode must be full, screening or download'
    def filters = file(params.filters, checkIfExists:true)
    if (!workflow.stubRun && params.mode != 'download') {
        ['carpdm_fasta','crass_fasta','reference_lock'].each { key ->
            if (!params[key]) error "Missing --${key}; see docs/INPUTS_NEEDED.md"
        }
    }
    if (params.mode == 'full') {
        ['kraken_db','amr_adapter','taxonomy_adapter','report_adapter'].each { key ->
            if (!params[key]) error "Full reproduction requires --${key}; see docs/INPUTS_NEEDED.md"
            file(params[key],checkIfExists:true)
        }
    }
    if (!(params.bowtie2_preset in ['very-sensitive','sensitive','very-sensitive-local','sensitive-local'])) {
        error 'Unsupported Bowtie2 preset'
    }
    def ena = params.ena_report ? file(params.ena_report,checkIfExists:true) : file("${projectDir}/assets/NO_ENA_REPORT")
    def refs = Channel.empty()
    def readyFilters = Channel.value(filters)
    if (params.mode != 'download') {
        def carp = file(params.carpdm_fasta ?: "${projectDir}/tests/reference.fa",checkIfExists:true)
        def crass = file(params.crass_fasta ?: "${projectDir}/tests/reference.fa",checkIfExists:true)
        def lock = file(params.reference_lock ?: "${projectDir}/assets/reference_lock.example.json",checkIfExists:true)
        VERIFY_REFERENCES(lock,carp,crass)
        readyFilters = VERIFY_REFERENCES.out.map { ignored -> filters }
        // Index building depends on verification; failures cannot launch mappings.
        refs = VERIFY_REFERENCES.out.map { ignored -> [['carpdm',carp],['crass',crass]] }.flatMap()
        BUILD_INDEX(refs)
    }
    RESOLVE_ENA(readyFilters,ena)
    downloads = RESOLVE_ENA.out.downloads.splitCsv(header:true,sep:'\t').map { r ->
        tuple(r.filter_id,r.run_accession,r.url1,r.url2,r.md5_1,r.md5_2)
    }
    DOWNLOAD(downloads)
    MERGE_READS(DOWNLOAD.out.groupTuple())
    if (params.mode == 'download') {
        MERGE_READS.out.view { id,reads -> "Verified downloaded reads: ${id}: ${reads}" }
    } else {
        FASTQC(MERGE_READS.out)
        FASTP(MERGE_READS.out)
        MULTIQC(FASTQC.out.mix(FASTP.out.json,FASTP.out.html).collect())
        MAP_SCREEN(FASTP.out.reads.combine(BUILD_INDEX.out))
        POOL_SCREEN(MAP_SCREEN.out.metrics.groupTuple(),FASTP.out.json.collect(),filters)
        if (params.mode == 'full') {
            AMR_ORIGINAL(FASTP.out.reads,file(params.amr_adapter))
            KRAKEN(FASTP.out.reads,file(params.kraken_db))
            TAXONOMY_ORIGINAL(KRAKEN.out,file(params.kraken_db),file(params.taxonomy_adapter))
            REPORT_ORIGINAL(file(params.report_adapter),filters,
                AMR_ORIGINAL.out.map { id,dir -> dir }.collect(),
                TAXONOMY_ORIGINAL.out.map { id,dir -> dir }.collect(),
                MAP_SCREEN.out.bam.map { id,label,bam -> bam }.collect(),
                MAP_SCREEN.out.metrics.map { label,metrics -> metrics }.collect(),
                POOL_SCREEN.out.map { label,path -> path }.collect(),
                FASTP.out.json.collect(),VERIFY_REFERENCES.out)
        }
    }
}
