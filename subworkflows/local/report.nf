/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT
    Comprehensive reporting with all QC metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FILTER_INFO_BY_TARGET     } from '../../modules/local/report/filter_info_by_target'
include { REPORT_WELL_IMPUTED       } from '../../modules/local/report/report_well_imputed'
include { PLOT_PERFORMANCE          } from '../../modules/local/report/plot_performance'
include { GENERATE_CHUNK_JSON       } from '../../modules/local/report/generate_chunk_json'
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
// New imputation analysis modules - RE-ENABLED AFTER FIXING SYNTAX ERRORS
include { COMPARE_PRE_POST_IMPUTATION } from '../../modules/local/report/compare_pre_post_imputation'
include { PLOT_R2_GENOMIC_WINDOWS     } from '../../modules/local/report/plot_r2_genomic_windows'
// New advanced QC modules - PARTIALLY ENABLED (heterozygosity and HWE have syntax issues)
include { PLOT_DOSAGE_DISTRIBUTION  } from '../../modules/local/report/plot_dosage_distribution'
include { PLOT_CALIBRATION          } from '../../modules/local/report/plot_calibration'
include { PLOT_CONCORDANCE_MAF      } from '../../modules/local/report/plot_concordance_maf'
include { PLOT_CROSS_VALIDATION     } from '../../modules/local/report/plot_cross_validation'
// include { PLOT_HETEROZYGOSITY       } from '../../modules/local/report/plot_heterozygosity'  // Syntax error - needs fixing
// include { PLOT_HWE_DEVIATION        } from '../../modules/local/report/plot_hwe_deviation'    // Syntax error - needs fixing
include { GENERATE_SUMMARY_REPORT   } from '../../modules/local/report/generate_summary'
// Terra-style visualization modules
include { PLOT_IMPUTATION_ACCURACY_MAF_BINS } from '../../modules/local/report/plot_imputation_accuracy_maf_bins'
include { PLOT_AGGREGATED_R2_DASHBOARD      } from '../../modules/local/report/plot_aggregated_r2_dashboard'
include { CONVERT_VCF_TO_INFO                } from '../../modules/local/report/convert_vcf_to_info'

workflow REPORT {
    take:
    ch_imputed // channel: [ val(meta), val(ref_name), path(vcf), path(vcf_index), path(info) ]
    ch_ref_vcf // channel: [ val(meta), val(ref_name), path(ref_vcf) ]
    ch_pre_imputation_vcf // channel: [ val(meta), path(vcf), path(vcf_index) ] - Original input VCF for comparison

    main:
    ch_versions = Channel.empty()
    
    //
    // MODULE: Convert VCF sites files to info format
    // This preprocessing step ensures all plotting modules receive .info files
    // MODIFIED: Add error tolerance for missing chunks
    //
    ch_vcf_sites = ch_imputed
        .map { meta, ref_name, vcf, vcf_index, info ->
            // Extract the sites VCF file for conversion
            // The 'info' here is actually the .sites.vcf.gz file from minimac4
            [meta, ref_name, info]
        }
    
    CONVERT_VCF_TO_INFO ( ch_vcf_sites )
    ch_versions = ch_versions.mix(CONVERT_VCF_TO_INFO.out.versions)
    
    // Now use the converted .info files for all downstream processing
    ch_info_only = CONVERT_VCF_TO_INFO.out.info
    
    //
    // MODULE: Filter info by target
    //
    FILTER_INFO_BY_TARGET ( ch_info_only )
    ch_versions = ch_versions.mix(FILTER_INFO_BY_TARGET.out.versions)
    ch_filtered_info = FILTER_INFO_BY_TARGET.out.filtered
    
    //
    // MODULE: Report well imputed variants
    // MODIFIED: Use error-tolerant channel
    //
    REPORT_WELL_IMPUTED ( ch_filtered_info )
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
    // MODIFIED: Use error-tolerant channel
    //
    REPORT_ACCURACY ( ch_filtered_info )
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
    // MODULE: Generate JSON summary for chunk aggregation
    //
    ch_chunk_json_input = REPORT_ACCURACY.out.report
        .join(REPORT_WELL_IMPUTED.out.report, by: [0, 1])
        .map { meta, ref_name, accuracy_txt, accuracy_tsv, well_imputed_txt, well_imputed_summary ->
            // Select the files we need for JSON generation
            [meta, ref_name, accuracy_txt, well_imputed_txt, well_imputed_summary]
        }
    
    GENERATE_CHUNK_JSON ( ch_chunk_json_input )
    ch_versions = ch_versions.mix(GENERATE_CHUNK_JSON.out.versions)
    
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
    // MODULE: Pre/Post-Imputation Comparison
    // Compare original input VCF with imputed output
    //
    ch_pre_post_comparison = ch_pre_imputation_vcf
        .join(ch_imputed.map { meta, ref_name, vcf, vcf_index, info -> 
            [meta, vcf, vcf_index] 
        })
        .map { meta, pre_vcf, pre_index, post_vcf, post_index ->
            [meta, pre_vcf, pre_index, post_vcf, post_index]
        }
    
    COMPARE_PRE_POST_IMPUTATION ( ch_pre_post_comparison )
    ch_versions = ch_versions.mix(COMPARE_PRE_POST_IMPUTATION.out.versions)
    
    //
    // MODULE: R2 Genomic Windows Analysis
    // Identify poorly imputed regions using sliding windows
    //
    PLOT_R2_GENOMIC_WINDOWS ( ch_info_only )
    ch_versions = ch_versions.mix(PLOT_R2_GENOMIC_WINDOWS.out.versions)
    
    //
    // MODULE: Advanced QC Plotting Modules
    //
    // Dosage distribution analysis
    //
    ch_dosage_input = ch_imputed
        .map { meta, ref_name, vcf, vcf_index, info ->
            [meta, ref_name, vcf, vcf_index]
        }
    
    PLOT_DOSAGE_DISTRIBUTION ( ch_dosage_input )
    ch_versions = ch_versions.mix(PLOT_DOSAGE_DISTRIBUTION.out.versions)
    
    // Calibration plot - DISABLED (requires validation/truth data)
    // NOTE: Requires external validation dataset to compare against
    // PLOT_CALIBRATION ( ch_info_only )
    // ch_versions = ch_versions.mix(PLOT_CALIBRATION.out.versions)
    
    // Concordance MAF plot - DISABLED (requires validation/truth data)
    // NOTE: Requires external validation dataset for concordance analysis
    // PLOT_CONCORDANCE_MAF ( ch_info_only )
    // ch_versions = ch_versions.mix(PLOT_CONCORDANCE_MAF.out.versions)
    
    // Cross-validation plot - DISABLED (requires validation/truth data)
    // NOTE: Requires held-out validation data from cross-validation run
    // PLOT_CROSS_VALIDATION ( ch_info_only )
    // ch_versions = ch_versions.mix(PLOT_CROSS_VALIDATION.out.versions)
    
    // Heterozygosity analysis - DISABLED due to syntax errors
    // ch_het_input = ch_imputed.map { meta, ref_name, vcf, vcf_index, info ->
    //     [meta, ref_name, vcf, vcf_index]
    // }
    // PLOT_HETEROZYGOSITY ( ch_het_input )
    // ch_versions = ch_versions.mix(PLOT_HETEROZYGOSITY.out.versions)
    
    // HWE deviation plot - DISABLED due to syntax errors  
    // ch_hwe_input = ch_imputed.map { meta, ref_name, vcf, vcf_index, info ->
    //     [meta, ref_name, vcf, vcf_index]
    // }
    // PLOT_HWE_DEVIATION ( ch_hwe_input )
    // ch_versions = ch_versions.mix(PLOT_HWE_DEVIATION.out.versions)
    
    // Terra-style MAF bins accuracy analysis
    PLOT_IMPUTATION_ACCURACY_MAF_BINS ( ch_info_only )
    ch_versions = ch_versions.mix(PLOT_IMPUTATION_ACCURACY_MAF_BINS.out.versions)
    
    // Aggregated R² dashboard - Now enabled with converted .info files
    PLOT_AGGREGATED_R2_DASHBOARD ( ch_info_only )
    ch_versions = ch_versions.mix(PLOT_AGGREGATED_R2_DASHBOARD.out.versions)
    
    // Generate summary report - DISABLED (needs input format fix)
    // ch_summary_input = REPORT_WELL_IMPUTED.out.report
    //     .join(REPORT_ACCURACY.out.report, by: [0, 1])
    //     .join(AVERAGE_R2.out.average, by: [0, 1])
    //     .map { meta, ref_name, well_txt, well_summary, acc_txt, acc_tsv, avg_r2 ->
    //         [meta, ref_name, well_summary, acc_txt, avg_r2]
    //     }
    // GENERATE_SUMMARY_REPORT ( ch_summary_input )
    // ch_versions = ch_versions.mix(GENERATE_SUMMARY_REPORT.out.versions)
    
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
                   AVERAGE_R2.out.summary,
                   GENERATE_CHUNK_JSON.out.json,
                   // GENERATE_SUMMARY_REPORT.out.report,  // Disabled - needs input format fix
                   COMPARE_PRE_POST_IMPUTATION.out.comparison
                   // MERGE_FREQ_COMPARISON.out.summary  // Disabled - outputs path not tuple, breaks mix
               )
    chunk_json = GENERATE_CHUNK_JSON.out.json  // Separate emit for aggregation
    well_imputed = REPORT_WELL_IMPUTED.out.report  // For backward compatibility
    plots    = ch_performance_plots.mix(
                   ch_accuracy_plots,
                   ch_r2_maf_plots,
                   ch_freq_comparison_plots,
                   ch_r2_snppos_plots,
                   ch_r2_snpcount_plots,
                   ch_hist_r2_snpcount_plots,
                   ch_maf_r2_plots,
                   PLOT_R2_GENOMIC_WINDOWS.out.plot,
                   PLOT_DOSAGE_DISTRIBUTION.out.plot,
                   // PLOT_CALIBRATION.out.plot,
                   // PLOT_CONCORDANCE_MAF.out.plot,
                   // PLOT_CROSS_VALIDATION.out.plot,
                   PLOT_IMPUTATION_ACCURACY_MAF_BINS.out.plot,
                   PLOT_AGGREGATED_R2_DASHBOARD.out.plot
                   // PLOT_HETEROZYGOSITY.out.plot,  // Disabled - syntax errors
                   // PLOT_HWE_DEVIATION.out.plot    // Disabled - syntax errors
               )
    // Specific outputs for hierarchical aggregation
    well_imputed = REPORT_WELL_IMPUTED.out.report
    performance_plots = ch_performance_plots
    filtered_info = FILTER_INFO_BY_TARGET.out.filtered
    versions = ch_versions
}