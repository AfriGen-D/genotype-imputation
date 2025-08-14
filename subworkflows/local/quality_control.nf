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
    
    // Remove duplicate variants
    REMOVE_DUPLICATES(ch_vcfs)
    ch_versions = ch_versions.mix(REMOVE_DUPLICATES.out.versions)
    
    // Split multi-allelic variants
    SPLIT_MULTIALLELIC(REMOVE_DUPLICATES.out.vcf)
    ch_versions = ch_versions.mix(SPLIT_MULTIALLELIC.out.versions)
    
    // Filter variants based on quality metrics
    FILTER_VARIANTS(
        SPLIT_MULTIALLELIC.out.vcf,
        qc_params.min_ac ?: 2,
        qc_params.min_maf ?: 0.01,
        qc_params.max_missing ?: 0.05
    )
    ch_versions = ch_versions.mix(FILTER_VARIANTS.out.versions)
    
    // Calculate sample and variant missingness
    CALCULATE_MISSINGNESS(FILTER_VARIANTS.out.vcf)
    ch_versions = ch_versions.mix(CALCULATE_MISSINGNESS.out.versions)
    
    // Hardy-Weinberg equilibrium check
    CHECK_HWE(
        FILTER_VARIANTS.out.vcf,
        qc_params.hwe_pvalue ?: 1e-6
    )
    ch_versions = ch_versions.mix(CHECK_HWE.out.versions)
    
    // Annotate variants with AF and other tags
    // Pass null as annotation_db if not available
    ANNOTATE_VARIANTS(
        FILTER_VARIANTS.out.vcf,
        file("NO_FILE")
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)
    
    // Collect QC metrics
    ch_qc_metrics = REMOVE_DUPLICATES.out.stats
        .join(SPLIT_MULTIALLELIC.out.stats)
        .join(FILTER_VARIANTS.out.stats)
        .join(CALCULATE_MISSINGNESS.out.report)
        .join(CHECK_HWE.out.stats)
        .join(ANNOTATE_VARIANTS.out.stats)
    
    emit:
    vcf      = ANNOTATE_VARIANTS.out.vcf
    metrics  = ch_qc_metrics
    versions = ch_versions
}