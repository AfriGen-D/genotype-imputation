/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FILTER_INFO_BY_TARGET     } from '../../modules/local/report/filter_info_by_target'
include { REPORT_WELL_IMPUTED       } from '../../modules/local/report/report_well_imputed'
include { PLOT_PERFORMANCE          } from '../../modules/local/report/plot_performance'
include { REPORT_ACCURACY           } from '../../modules/local/report/report_accuracy'
include { PLOT_ACCURACY             } from '../../modules/local/report/plot_accuracy'
include { PLOT_R2_MAF               } from '../../modules/local/report/plot_r2_maf'
include { PLOT_FREQ_COMPARISON      } from '../../modules/local/report/plot_freq_comparison'
include { COLLECT_WARNINGS          } from '../../modules/local/report/collect_warnings'

workflow REPORT {
    take:
    ch_imputed // channel: [ val(meta), path(vcf), path(info) ]

    main:
    ch_versions = Channel.empty()
    
    //
    // MODULE: Filter info by target
    //
    FILTER_INFO_BY_TARGET ( ch_imputed )
    ch_versions = ch_versions.mix(FILTER_INFO_BY_TARGET.out.versions)
    
    //
    // MODULE: Report well imputed variants
    //
    REPORT_WELL_IMPUTED ( FILTER_INFO_BY_TARGET.out.filtered )
    ch_versions = ch_versions.mix(REPORT_WELL_IMPUTED.out.versions)
    
    //
    // MODULE: Plot imputation performance
    //
    PLOT_PERFORMANCE ( REPORT_WELL_IMPUTED.out.report )
    ch_versions = ch_versions.mix(PLOT_PERFORMANCE.out.versions)
    
    //
    // MODULE: Report accuracy metrics
    //
    REPORT_ACCURACY ( FILTER_INFO_BY_TARGET.out.filtered )
    ch_versions = ch_versions.mix(REPORT_ACCURACY.out.versions)
    
    //
    // MODULE: Plot accuracy
    //
    PLOT_ACCURACY ( REPORT_ACCURACY.out.report )
    ch_versions = ch_versions.mix(PLOT_ACCURACY.out.versions)
    
    //
    // MODULE: Plot R2 vs MAF
    //
    PLOT_R2_MAF ( FILTER_INFO_BY_TARGET.out.filtered )
    ch_versions = ch_versions.mix(PLOT_R2_MAF.out.versions)
    
    //
    // MODULE: Plot frequency comparison
    //
    PLOT_FREQ_COMPARISON ( ch_imputed )
    ch_versions = ch_versions.mix(PLOT_FREQ_COMPARISON.out.versions)
    
    //
    // MODULE: Collect warnings from pipeline log
    //
    // Create a channel with the Nextflow log file
    ch_log = ch_imputed
        .first()
        .map { meta, ref_name, info ->
            [meta, file("${workflow.launchDir}/.nextflow.log")]
        }
    
    COLLECT_WARNINGS ( ch_log )
    ch_versions = ch_versions.mix(COLLECT_WARNINGS.out.versions)
    
    emit:
    reports  = REPORT_WELL_IMPUTED.out.report.mix(
                   REPORT_ACCURACY.out.report,
                   COLLECT_WARNINGS.out.warnings,
                   COLLECT_WARNINGS.out.summary
               )
    plots    = PLOT_PERFORMANCE.out.plot.mix(
                   PLOT_ACCURACY.out.plot,
                   PLOT_R2_MAF.out.plot,
                   PLOT_FREQ_COMPARISON.out.plot
               )
    versions = ch_versions
}