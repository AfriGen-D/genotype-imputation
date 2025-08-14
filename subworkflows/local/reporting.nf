/*
========================================================================================
    SUBWORKFLOW: REPORTING
========================================================================================
*/

include { CALCULATE_IMPUTATION_METRICS } from '../../modules/local/visualization/calculate_metrics'
include { GENERATE_PERFORMANCE_REPORT } from '../../modules/local/visualization/performance_report'
include { GENERATE_ACCURACY_REPORT } from '../../modules/local/visualization/accuracy_report'
include { PLOT_MAF_CONCORDANCE } from '../../modules/local/visualization/plot_maf_concordance'
include { PLOT_INFO_DISTRIBUTION } from '../../modules/local/visualization/plot_info_distribution'
include { PLOT_R2_BY_POSITION } from '../../modules/local/visualization/plot_r2_position'
include { CREATE_HTML_REPORT } from '../../modules/local/visualization/html_report'
include { GENERATE_MULTIQC_CONFIG } from '../../modules/local/visualization/multiqc_config'

workflow REPORTING {
    take:
    ch_imputed_vcfs   // channel: [ val(meta), path(vcf) ]
    ch_info_scores    // channel: [ val(meta), path(info) ]
    ch_qc_metrics     // channel: [ val(meta), path(metrics) ]
    ch_original_vcfs  // channel: [ val(meta), path(vcf) ] - for concordance
    report_level      // val: 'summary', 'detailed', or 'full'

    main:
    ch_versions = Channel.empty()
    
    // Calculate comprehensive imputation metrics
    ch_vcf_with_info = ch_imputed_vcfs
        .join(ch_info_scores)
    
    CALCULATE_IMPUTATION_METRICS(
        ch_vcf_with_info,
        ch_original_vcfs
    )
    ch_versions = ch_versions.mix(CALCULATE_IMPUTATION_METRICS.out.versions)
    
    // Generate performance report
    GENERATE_PERFORMANCE_REPORT(
        CALCULATE_IMPUTATION_METRICS.out.metrics,
        report_level
    )
    ch_versions = ch_versions.mix(GENERATE_PERFORMANCE_REPORT.out.versions)
    
    // Generate accuracy report by MAF bins
    GENERATE_ACCURACY_REPORT(
        CALCULATE_IMPUTATION_METRICS.out.metrics,
        ch_info_scores
    )
    ch_versions = ch_versions.mix(GENERATE_ACCURACY_REPORT.out.versions)
    
    // Create visualizations
    // Create empty channel for optional masked file input
    ch_empty_masked = Channel.empty()
    
    PLOT_MAF_CONCORDANCE(
        CALCULATE_IMPUTATION_METRICS.out.concordance_data,
        ch_empty_masked
    )
    ch_versions = ch_versions.mix(PLOT_MAF_CONCORDANCE.out.versions)
    
    PLOT_INFO_DISTRIBUTION(
        ch_info_scores
    )
    ch_versions = ch_versions.mix(PLOT_INFO_DISTRIBUTION.out.versions)
    
    PLOT_R2_BY_POSITION(
        CALCULATE_IMPUTATION_METRICS.out.r2_data,
        ch_empty_masked
    )
    ch_versions = ch_versions.mix(PLOT_R2_BY_POSITION.out.versions)
    
    // Combine all reports and plots
    ch_all_reports = GENERATE_PERFORMANCE_REPORT.out.report
        .join(GENERATE_ACCURACY_REPORT.out.report)
        .join(PLOT_MAF_CONCORDANCE.out.plot)
        .join(PLOT_INFO_DISTRIBUTION.out.plot)
        .join(PLOT_R2_BY_POSITION.out.plot)
    
    // Generate comprehensive HTML report
    CREATE_HTML_REPORT(
        CALCULATE_IMPUTATION_METRICS.out.metrics,
        ch_all_reports,
        ch_qc_metrics,
        "Imputation Quality Report - ${report_level}"
    )
    ch_versions = ch_versions.mix(CREATE_HTML_REPORT.out.versions)
    
    // Generate MultiQC configuration
    if (report_level == 'full') {
        GENERATE_MULTIQC_CONFIG(
            CREATE_HTML_REPORT.out.html_report.collect(),
            ch_versions.unique().collect()
        )
        ch_final_report = GENERATE_MULTIQC_CONFIG.out.report
    } else {
        ch_final_report = CREATE_HTML_REPORT.out.html_report
    }
    
    emit:
    html_report  = ch_final_report
    plots        = PLOT_MAF_CONCORDANCE.out.plot.mix(
                     PLOT_INFO_DISTRIBUTION.out.plot,
                     PLOT_R2_BY_POSITION.out.plot
                   )
    metrics      = CALCULATE_IMPUTATION_METRICS.out.metrics
    versions     = ch_versions
}