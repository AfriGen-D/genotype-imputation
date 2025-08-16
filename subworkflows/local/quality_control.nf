/*
========================================================================================
    SUBWORKFLOW: QUALITY_CONTROL
========================================================================================
*/

include { REMOVE_DUPLICATES } from '../../modules/local/qc/remove_duplicates'
include { SPLIT_MULTIALLELIC } from '../../modules/local/qc/split_multiallelic'
include { FILTER_VARIANTS } from '../../modules/local/qc/filter_variants'
include { CALCULATE_MISSINGNESS } from '../../modules/local/qc/calculate_missingness'
include { CHECK_HWE } from '../../modules/local/qc/check_hwe'
include { ANNOTATE_VARIANTS } from '../../modules/local/qc/annotate_variants'

workflow QUALITY_CONTROL {
    take:
    ch_vcfs        // channel: [ val(meta), path(vcf) ]
    qc_params      // map: QC parameters

    main:
    ch_versions = Channel.empty()
    ch_current_vcf = ch_vcfs
    ch_qc_metrics = Channel.empty()
    
    // ESSENTIAL: Remove duplicate variants (always run)
    if (params.remove_duplicates) {
        REMOVE_DUPLICATES(ch_current_vcf)
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES.out.versions)
        ch_current_vcf = REMOVE_DUPLICATES.out.vcf
        ch_qc_metrics = ch_qc_metrics.mix(REMOVE_DUPLICATES.out.stats)
    }
    
    // ESSENTIAL: Split multi-allelic variants (always run for correct imputation)
    if (params.split_multiallelic) {
        SPLIT_MULTIALLELIC(ch_current_vcf)
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELIC.out.versions)
        ch_current_vcf = SPLIT_MULTIALLELIC.out.vcf
        ch_qc_metrics = ch_qc_metrics.mix(SPLIT_MULTIALLELIC.out.stats)
    }
    
    // ESSENTIAL: Filter variants based on quality metrics
    FILTER_VARIANTS(
        ch_current_vcf,
        qc_params.min_ac ?: 2,
        qc_params.min_maf ?: 0.01,
        qc_params.max_missing ?: 0.05
    )
    ch_versions = ch_versions.mix(FILTER_VARIANTS.out.versions)
    ch_current_vcf = FILTER_VARIANTS.out.vcf
    ch_qc_metrics = ch_qc_metrics.mix(FILTER_VARIANTS.out.stats)
    
    // OPTIONAL: Calculate missingness (skip in minimal mode)
    if (!params.skip_missingness_calc && !params.minimal_qc_mode) {
        CALCULATE_MISSINGNESS(ch_current_vcf)
        ch_versions = ch_versions.mix(CALCULATE_MISSINGNESS.out.versions)
        ch_qc_metrics = ch_qc_metrics.mix(CALCULATE_MISSINGNESS.out.report)
    }
    
    // OPTIONAL: Hardy-Weinberg equilibrium check (skip in minimal mode)
    if (!params.skip_hwe_check && !params.minimal_qc_mode) {
        CHECK_HWE(
            ch_current_vcf,
            qc_params.hwe_pvalue ?: 1e-6
        )
        ch_versions = ch_versions.mix(CHECK_HWE.out.versions)
        ch_qc_metrics = ch_qc_metrics.mix(CHECK_HWE.out.stats)
    }
    
    // OPTIONAL: Annotate variants (skip in minimal mode)
    if (!params.skip_annotation && !params.minimal_qc_mode) {
        ANNOTATE_VARIANTS(
            ch_current_vcf,
            file("NO_FILE")
        )
        ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)
        ch_current_vcf = ANNOTATE_VARIANTS.out.vcf
        ch_qc_metrics = ch_qc_metrics.mix(ANNOTATE_VARIANTS.out.stats)
    }
    
    emit:
    vcf      = ch_current_vcf
    metrics  = ch_qc_metrics
    versions = ch_versions
}