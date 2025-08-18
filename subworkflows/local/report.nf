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
include { EXTRACT_ALLELE_FREQ       } from '../../modules/local/report/extract_allele_freq'
include { COMBINE_FREQ_BY_CHR       } from '../../modules/local/report/combine_freq_by_chr'
include { COMBINE_FREQ_GENOME       } from '../../modules/local/report/combine_freq_genome'
include { PLOT_FREQ_COMPARISON      } from '../../modules/local/report/plot_freq_comparison'
include { COLLECT_WARNINGS          } from '../../modules/local/report/collect_warnings'
include { PLOT_R2_SNPPOS            } from '../../modules/local/report/plot_r2_snppos'
include { PLOT_R2_SNPCOUNT          } from '../../modules/local/report/plot_r2_snpcount'
include { PLOT_HIST_R2_SNPCOUNT     } from '../../modules/local/report/plot_hist_r2_snpcount'
include { PLOT_MAF_R2               } from '../../modules/local/report/plot_maf_r2'
include { AVERAGE_R2                } from '../../modules/local/report/average_r2'
// Import frequency extraction modules
include { EXTRACT_FREQ_SIMPLE       } from '../../modules/local/report/extract_freq_simple'
include { EXTRACT_FREQ_COMPARISON   } from '../../modules/local/report/extract_freq_comparison'
include { EXTRACT_FREQ_IMPUTED      } from '../../modules/local/report/extract_freq_imputed'
include { EXTRACT_FREQ_REFERENCE    } from '../../modules/local/report/extract_freq_reference'
include { MERGE_FREQ_COMPARISON     } from '../../modules/local/report/merge_freq_comparison'
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
    ch_imputed // channel: [ val(meta), val(ref_name), path(vcf), path(vcf_index), path(info) ]
    ch_ref_vcf // channel: [ val(meta), val(ref_name), path(ref_vcf) ]

    main:
    ch_versions = Channel.empty()
    
    //
    // MODULE: Filter info by target
    //
    // Extract just the info file for FILTER_INFO_BY_TARGET
    ch_info_only = ch_imputed.map { meta, ref_name, vcf, vcf_index, info ->
        [meta, ref_name, info]
    }
    FILTER_INFO_BY_TARGET ( ch_info_only )
    ch_versions = ch_versions.mix(FILTER_INFO_BY_TARGET.out.versions)
    
    //
    // MODULE: Report well imputed variants
    //
    REPORT_WELL_IMPUTED ( FILTER_INFO_BY_TARGET.out.filtered )
    ch_versions = ch_versions.mix(REPORT_WELL_IMPUTED.out.versions)
    
    //
    // MODULE: Plot imputation performance (chunk-level)
    // Conditional execution based on params.generate_chunk_plots
    //
    if (params.generate_chunk_plots) {
        PLOT_PERFORMANCE ( REPORT_WELL_IMPUTED.out.report )
        ch_versions = ch_versions.mix(PLOT_PERFORMANCE.out.versions)
        ch_performance_plots = PLOT_PERFORMANCE.out.plot
    } else {
        ch_performance_plots = Channel.empty()
    }
    
    //
    // MODULE: Report accuracy metrics
    //
    REPORT_ACCURACY ( FILTER_INFO_BY_TARGET.out.filtered )
    ch_versions = ch_versions.mix(REPORT_ACCURACY.out.versions)
    
    //
    // MODULE: Plot accuracy (chunk-level)
    // Conditional execution based on params.generate_chunk_plots
    //
    if (params.generate_chunk_plots) {
        PLOT_ACCURACY ( REPORT_ACCURACY.out.report )
        ch_versions = ch_versions.mix(PLOT_ACCURACY.out.versions)
        ch_accuracy_plots = PLOT_ACCURACY.out.plot
    } else {
        ch_accuracy_plots = Channel.empty()
    }
    
    //
    // MODULE: Plot R2 vs MAF (chunk-level)
    // Conditional execution based on params.generate_chunk_plots
    //
    if (params.generate_chunk_plots) {
        PLOT_R2_MAF ( FILTER_INFO_BY_TARGET.out.filtered )
        ch_versions = ch_versions.mix(PLOT_R2_MAF.out.versions)
        ch_r2_maf_plots = PLOT_R2_MAF.out.plot
    } else {
        ch_r2_maf_plots = Channel.empty()
    }
    
    //
    // MODULE: Frequency comparison pipeline (modular approach)
    //
    // Step 1: Extract frequencies from imputed VCFs
    ch_imputed_for_freq = ch_imputed.map { meta, ref_name, vcf, vcf_index, info ->
        [meta, ref_name, vcf, vcf_index]
    }
    
    EXTRACT_FREQ_IMPUTED ( ch_imputed_for_freq )
    ch_versions = ch_versions.mix(EXTRACT_FREQ_IMPUTED.out.versions)
    
    // Step 2: Extract frequencies from reference panel VCFs (if available)
    // Use left join to handle missing reference VCFs
    ch_ref_for_freq = ch_ref_vcf.map { meta, ref_name, ref_vcf ->
        [meta, ref_name, ref_vcf ?: file("NO_FILE")]
    }
    
    EXTRACT_FREQ_REFERENCE ( ch_ref_for_freq )
    ch_versions = ch_versions.mix(EXTRACT_FREQ_REFERENCE.out.versions)
    
    // Step 3: Merge and compare frequencies
    ch_freq_to_merge = EXTRACT_FREQ_IMPUTED.out.frequencies
        .join(EXTRACT_FREQ_REFERENCE.out.frequencies, by: [0, 1], remainder: true)
        .map { meta, ref_name, imp_freq, ref_freq ->
            // Handle missing reference frequencies
            def ref_freq_file = ref_freq ?: file("NO_FILE")
            [meta, ref_name, imp_freq, ref_freq_file]
        }
    
    MERGE_FREQ_COMPARISON ( ch_freq_to_merge )
    ch_versions = ch_versions.mix(MERGE_FREQ_COMPARISON.out.versions)
    
    // Step 4: Create frequency comparison plots (chunk-level)
    // Conditional execution based on params.generate_chunk_plots
    ch_freq_for_plot = MERGE_FREQ_COMPARISON.out.comparison
    
    if (params.generate_chunk_plots) {
        PLOT_FREQ_COMPARISON ( ch_freq_for_plot )
        ch_versions = ch_versions.mix(PLOT_FREQ_COMPARISON.out.versions)
        ch_freq_comparison_plots = PLOT_FREQ_COMPARISON.out.plot
    } else {
        ch_freq_comparison_plots = Channel.empty()
    }
    
    //
    // MODULE: Additional R² analysis plots (chunk-level)
    // Conditional execution based on params.generate_chunk_plots
    //
    if (params.generate_chunk_plots) {
        PLOT_R2_SNPPOS ( ch_info_only )
        ch_versions = ch_versions.mix(PLOT_R2_SNPPOS.out.versions)
        ch_r2_snppos_plots = PLOT_R2_SNPPOS.out.plot
        
        PLOT_R2_SNPCOUNT ( ch_info_only )
        ch_versions = ch_versions.mix(PLOT_R2_SNPCOUNT.out.versions)
        ch_r2_snpcount_plots = PLOT_R2_SNPCOUNT.out.plot
        
        PLOT_HIST_R2_SNPCOUNT ( ch_info_only )
        ch_versions = ch_versions.mix(PLOT_HIST_R2_SNPCOUNT.out.versions)
        ch_hist_r2_snpcount_plots = PLOT_HIST_R2_SNPCOUNT.out.plot
        
        PLOT_MAF_R2 ( ch_info_only )
        ch_versions = ch_versions.mix(PLOT_MAF_R2.out.versions)
        ch_maf_r2_plots = PLOT_MAF_R2.out.plot
    } else {
        ch_r2_snppos_plots = Channel.empty()
        ch_r2_snpcount_plots = Channel.empty()
        ch_hist_r2_snpcount_plots = Channel.empty()
        ch_maf_r2_plots = Channel.empty()
    }
    
    //
    // MODULE: Calculate average R²
    //
    AVERAGE_R2 ( ch_info_only )
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
                   // MERGE_FREQ_COMPARISON.out.summary  // Disabled - outputs path not tuple, breaks mix
               )
    plots    = ch_performance_plots.mix(
                   ch_accuracy_plots,
                   ch_r2_maf_plots,
                   ch_freq_comparison_plots,
                   ch_r2_snppos_plots,
                   ch_r2_snpcount_plots,
                   ch_hist_r2_snpcount_plots,
                   ch_maf_r2_plots
               )
    // Specific outputs for hierarchical aggregation
    well_imputed = REPORT_WELL_IMPUTED.out.report
    performance_plots = ch_performance_plots
    filtered_info = FILTER_INFO_BY_TARGET.out.filtered
    versions = ch_versions
}