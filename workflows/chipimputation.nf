/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS            } from '../subworkflows/local/preprocess'
include { PHASE                 } from '../subworkflows/local/phase'
include { IMPUTE                } from '../subworkflows/local/impute'
include { REPORT                } from '../subworkflows/local/report'
// include { QC_PLOTS              } from '../subworkflows/local/qc_plots'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CHIPIMPUTATION {

    take:
    ch_input // channel: samplesheet read in from --input

    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // SUBWORKFLOW: Preprocessing and QC
    //
    PREPROCESS(
        ch_input
    )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)
    
    //
    // SUBWORKFLOW: Phasing
    //
    PHASE(
        PREPROCESS.out.vcf
    )
    ch_versions = ch_versions.mix(PHASE.out.versions)
    
    //
    // SUBWORKFLOW: Imputation
    //
    IMPUTE(
        PHASE.out.phased
    )
    ch_versions = ch_versions.mix(IMPUTE.out.versions)
    
    //
    // SUBWORKFLOW: Reporting
    //
    // Combine imputed VCFs with their info files for reporting
    ch_report_input = IMPUTE.out.imputed
        .join(IMPUTE.out.info, by: [0])  // Join by meta
        .map { meta, ref_name, vcf, vcf_index, info ->
            [meta, ref_name, info]  // Extract only what FILTER_INFO_BY_TARGET needs
        }
    
    REPORT(
        ch_report_input
    )
    ch_versions = ch_versions.mix(REPORT.out.versions)
    
    //
    // SUBWORKFLOW: QC Plots (disabled for now - needs implementation)
    //
    // if (params.qc_plots) {
    //     QC_PLOTS(
    //         IMPUTE.out.imputed
    //     )
    //     ch_versions = ch_versions.mix(QC_PLOTS.out.versions)
    //     ch_multiqc_files = ch_multiqc_files.mix(QC_PLOTS.out.plots)
    // }

    emit:
    imputed      = IMPUTE.out.imputed       // channel: [ val(meta), path(vcf) ]
    reports      = REPORT.out.reports       // channel: [ val(meta), path(reports) ]
    versions     = ch_versions              // channel: [ path(versions.yml) ]
    multiqc      = ch_multiqc_files        // channel: [ path(files) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/