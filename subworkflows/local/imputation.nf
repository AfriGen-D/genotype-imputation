/*
========================================================================================
    SUBWORKFLOW: IMPUTATION
========================================================================================
*/

include { MINIMAC4_IMPUTE } from '../../modules/local/imputation/minimac4'
include { IMPUTE5_IMPUTE } from '../../modules/local/imputation/impute5'
include { BEAGLE5_IMPUTE } from '../../modules/local/imputation/beagle5'
include { EXTRACT_INFO_SCORES } from '../../modules/local/imputation/extract_info'
include { FILTER_IMPUTED } from '../../modules/local/imputation/filter_imputed'
include { MERGE_IMPUTED_CHUNKS } from '../../modules/local/vcf/merge_imputed'

workflow IMPUTATION {
    take:
    ch_phased_vcfs    // channel: [ val(meta), path(vcf) ]
    ch_reference      // channel: [ val(ref_name), path(m3vcf/bcf), path(vcf) ]
    imputation_tool   // val: 'minimac4', 'impute5', or 'beagle5'
    impute_params     // map: imputation parameters

    main:
    ch_versions = Channel.empty()
    
    // Combine phased VCFs with reference panels
    ch_input_with_ref = ch_phased_vcfs
        .combine(ch_reference)
    
    // Perform imputation with selected tool
    if (imputation_tool == 'minimac4') {
        MINIMAC4_IMPUTE(
            ch_input_with_ref,
            impute_params.window ?: 500000,
            impute_params.min_ratio ?: 0.00001
        )
        ch_imputed = MINIMAC4_IMPUTE.out.imputed_vcf
        ch_info = MINIMAC4_IMPUTE.out.info
        ch_versions = ch_versions.mix(MINIMAC4_IMPUTE.out.versions)
        
    } else if (imputation_tool == 'impute5') {
        IMPUTE5_IMPUTE(
            ch_input_with_ref,
            impute_params.ne ?: 20000,
            impute_params.buffer ?: 250000
        )
        ch_imputed = IMPUTE5_IMPUTE.out.imputed_vcf
        ch_info = IMPUTE5_IMPUTE.out.info
        ch_versions = ch_versions.mix(IMPUTE5_IMPUTE.out.versions)
        
    } else if (imputation_tool == 'beagle5') {
        BEAGLE5_IMPUTE(
            ch_input_with_ref,
            impute_params.window ?: 40,
            impute_params.overlap ?: 4
        )
        ch_imputed = BEAGLE5_IMPUTE.out.imputed_vcf
        ch_info = BEAGLE5_IMPUTE.out.info
        ch_versions = ch_versions.mix(BEAGLE5_IMPUTE.out.versions)
        
    } else {
        error "Unknown imputation tool: ${imputation_tool}"
    }
    
    // Extract INFO scores from imputed VCFs
    EXTRACT_INFO_SCORES(
        ch_imputed,
        impute_params.info_cutoff ?: 0.3
    )
    ch_versions = ch_versions.mix(EXTRACT_INFO_SCORES.out.versions)
    
    // Filter imputed variants by INFO score
    ch_imputed_with_info = ch_imputed
        .join(EXTRACT_INFO_SCORES.out.info)
    
    FILTER_IMPUTED(
        ch_imputed,
        impute_params.info_cutoff ?: 0.3,
        0.01  // MAF cutoff
    )
    ch_versions = ch_versions.mix(FILTER_IMPUTED.out.versions)
    
    // Group filtered chunks by sample for merging
    ch_grouped_imputed = FILTER_IMPUTED.out.filtered_vcf
        .map { meta, vcf, ref_name -> 
            [ meta.id, meta, vcf ]
        }
        .groupTuple(by: 0)
        .map { id, metas, vcfs -> 
            [ metas[0], vcfs ]
        }
    
    // Merge imputed chunks if needed
    MERGE_IMPUTED_CHUNKS(ch_grouped_imputed)
    ch_versions = ch_versions.mix(MERGE_IMPUTED_CHUNKS.out.versions)
    
    emit:
    imputed_vcf = MERGE_IMPUTED_CHUNKS.out.merged_vcf
    info_scores = EXTRACT_INFO_SCORES.out.info
    metrics     = FILTER_IMPUTED.out.stats
    versions    = ch_versions
}