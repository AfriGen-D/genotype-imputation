/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT
    Comprehensive reporting with all QC metrics
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
include { PLOT_R2_SNPPOS            } from '../../modules/local/report/plot_r2_snppos'
include { PLOT_R2_SNPCOUNT          } from '../../modules/local/report/plot_r2_snpcount'
include { PLOT_HIST_R2_SNPCOUNT     } from '../../modules/local/report/plot_hist_r2_snpcount'
include { PLOT_MAF_R2               } from '../../modules/local/report/plot_maf_r2'
include { AVERAGE_R2                } from '../../modules/local/report/average_r2'
// New advanced QC modules
// include { PLOT_DOSAGE_DISTRIBUTION  } from '../../modules/local/report/plot_dosage_distribution'
// include { PLOT_CALIBRATION          } from '../../modules/local/report/plot_calibration'
// include { PLOT_CONCORDANCE_MAF      } from '../../modules/local/report/plot_concordance_maf'
// include { PLOT_CROSS_VALIDATION     } from '../../modules/local/report/plot_cross_validation'
// include { PLOT_HETEROZYGOSITY       } from '../../modules/local/report/plot_heterozygosity'
// include { PLOT_HWE_DEVIATION        } from '../../modules/local/report/plot_hwe_deviation'
// include { GENERATE_SUMMARY_REPORT   } from '../../modules/local/report/generate_summary'

workflow REPORT {
    take:
    ch_imputed // channel: [ val(meta), val(ref_name), path(info) ]

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
    // MODULE: Plot frequency comparison - Skip for now due to missing VCF
    //
    // PLOT_FREQ_COMPARISON requires the actual VCF file which we don't have in this channel
    // Would need to join with the imputed VCF channel from IMPUTE workflow
    // PLOT_FREQ_COMPARISON ( ch_imputed )
    // ch_versions = ch_versions.mix(PLOT_FREQ_COMPARISON.out.versions)
    
    //
    // MODULE: Additional R² analysis plots
    //
    PLOT_R2_SNPPOS ( ch_imputed )
    ch_versions = ch_versions.mix(PLOT_R2_SNPPOS.out.versions)
    
    PLOT_R2_SNPCOUNT ( ch_imputed )
    ch_versions = ch_versions.mix(PLOT_R2_SNPCOUNT.out.versions)
    
    PLOT_HIST_R2_SNPCOUNT ( ch_imputed )
    ch_versions = ch_versions.mix(PLOT_HIST_R2_SNPCOUNT.out.versions)
    
    PLOT_MAF_R2 ( ch_imputed )
    ch_versions = ch_versions.mix(PLOT_MAF_R2.out.versions)
    
    //
    // MODULE: Calculate average R²
    //
    AVERAGE_R2 ( ch_imputed )
    ch_versions = ch_versions.mix(AVERAGE_R2.out.versions)
    
    //
    // MODULE: Collect warnings from pipeline log
    //
    // Skip COLLECT_WARNINGS for now as it requires access to .nextflow.log
    // which may not be accessible during execution
    // TODO: Implement alternative warning collection mechanism
    // ch_log = ch_imputed
    //     .first()
    //     .map { meta, ref_name, info ->
    //         [meta, file("${workflow.launchDir}/.nextflow.log")]
    //     }
    // 
    // COLLECT_WARNINGS ( ch_log )
    // ch_versions = ch_versions.mix(COLLECT_WARNINGS.out.versions)
    
    emit:
    reports  = REPORT_WELL_IMPUTED.out.report.mix(
                   REPORT_ACCURACY.out.report,
                   // COLLECT_WARNINGS.out.warnings,  // Disabled - needs alternative implementation
                   // COLLECT_WARNINGS.out.summary,   // Disabled - needs alternative implementation
                   AVERAGE_R2.out.average,
                   AVERAGE_R2.out.summary
               )
    plots    = PLOT_PERFORMANCE.out.plot.mix(
                   PLOT_ACCURACY.out.plot,
                   PLOT_R2_MAF.out.plot,
                   // PLOT_FREQ_COMPARISON.out.plot,  // Disabled - needs VCF file
                   PLOT_R2_SNPPOS.out.plot,
                   PLOT_R2_SNPCOUNT.out.plot,
                   PLOT_HIST_R2_SNPCOUNT.out.plot,
                   PLOT_MAF_R2.out.plot
               )
    versions = ch_versions
}