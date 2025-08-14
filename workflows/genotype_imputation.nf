/*
========================================================================================
    IMPORT MODULES AND SUBWORKFLOWS
========================================================================================
*/

// Core subworkflows
include { INPUT_VALIDATION } from '../subworkflows/local/input_validation'
include { QUALITY_CONTROL } from '../subworkflows/local/quality_control'
include { PHASING } from '../subworkflows/local/phasing'
include { IMPUTATION } from '../subworkflows/local/imputation'
include { REPORTING } from '../subworkflows/local/reporting'

// Utility modules
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow GENOTYPE_IMPUTATION {

    // Check mandatory parameters
    if (!params.input) { 
        log.error "ERROR: Input samplesheet not specified with --input"
        System.exit(1)
    }
    
    if (!params.reference_panels || params.reference_panels.size() == 0) {
        log.error "ERROR: Reference panels not specified with --reference_panels"
        System.exit(1)
    }
    
    // Print workflow summary
    log.info """\
    ====================================
    GENOTYPE IMPUTATION PIPELINE v2.0
    ====================================
    Input            : ${params.input}
    Output directory : ${params.outdir}
    Genome build     : ${params.genome_build}
    Phasing tool     : ${params.phasing_tool}
    Imputation tool  : ${params.imputation_tool}
    Report level     : ${params.report_level}
    ====================================
    """.stripIndent()

    // Track software versions
    ch_versions = Channel.empty()
    
    /*
    ========================================================================================
        PREPARE INPUT CHANNELS
    ========================================================================================
    */
    
    // Create input channel from samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.population = row.population ?: 'Unknown'
            meta.sex = row.sex ?: 'Unknown'
            [meta, file(row.vcf)]
        }
        .set { ch_input_vcfs }
    
    // Create reference panel channel
    Channel
        .from(params.reference_panels)
        .map { ref ->
            def (name, m3vcf, vcf) = ref
            [name, file(m3vcf), file(vcf)]
        }
        .set { ch_reference_panels }
    
    /*
    ========================================================================================
        SUBWORKFLOW: INPUT VALIDATION
    ========================================================================================
    */
    
    INPUT_VALIDATION(
        ch_input_vcfs,
        ch_reference_panels,
        params.genome_build
    )
    ch_versions = ch_versions.mix(INPUT_VALIDATION.out.versions)
    
    /*
    ========================================================================================
        SUBWORKFLOW: QUALITY CONTROL
    ========================================================================================
    */
    
    def qc_params = [
        min_ac: params.qc_min_ac,
        min_maf: params.qc_min_maf,
        max_missing: params.qc_max_missing,
        hwe_pvalue: params.qc_hwe_pvalue
    ]
    
    QUALITY_CONTROL(
        INPUT_VALIDATION.out.vcf,
        qc_params
    )
    ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)
    
    /*
    ========================================================================================
        SUBWORKFLOW: PHASING
    ========================================================================================
    */
    
    PHASING(
        QUALITY_CONTROL.out.vcf,
        ch_reference_panels,
        params.eagle_genetic_map ? file(params.eagle_genetic_map) : null,
        params.phasing_tool,
        params.chunk_size
    )
    ch_versions = ch_versions.mix(PHASING.out.versions)
    
    /*
    ========================================================================================
        SUBWORKFLOW: IMPUTATION
    ========================================================================================
    */
    
    def impute_params = [
        window: params.impute_window,
        min_ratio: params.impute_min_ratio,
        ne: params.impute_ne,
        buffer: params.impute_buffer,
        info_cutoff: params.impute_info_cutoff
    ]
    
    IMPUTATION(
        PHASING.out.phased_vcf,
        ch_reference_panels,
        params.imputation_tool,
        impute_params
    )
    ch_versions = ch_versions.mix(IMPUTATION.out.versions)
    
    /*
    ========================================================================================
        SUBWORKFLOW: REPORTING
    ========================================================================================
    */
    
    // Prepare original VCFs for concordance analysis if requested
    ch_original_vcfs = params.concordance_analysis ? ch_input_vcfs : Channel.empty()
    
    REPORTING(
        IMPUTATION.out.imputed_vcf,
        IMPUTATION.out.info_scores,
        QUALITY_CONTROL.out.metrics,
        ch_original_vcfs,
        params.report_level
    )
    ch_versions = ch_versions.mix(REPORTING.out.versions)
    
    /*
    ========================================================================================
        MODULE: DUMP SOFTWARE VERSIONS
    ========================================================================================
    */
    
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
    emit:
    versions = ch_versions
    
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}