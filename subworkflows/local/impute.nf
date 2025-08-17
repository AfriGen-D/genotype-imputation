/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: IMPUTE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IMPUTE_MINIMAC4           } from '../../modules/local/impute/impute_minimac4'
// include { EXTRACT_IMPUTE_INFO       } from '../../modules/local/impute/extract_impute_info'
// include { COMBINE_IMPUTE            } from '../../modules/local/impute/combine_impute'
// include { COMBINE_INFO              } from '../../modules/local/impute/combine_info'

workflow IMPUTE {
    take:
    ch_phased // channel: [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()
    
    // Simplified imputation - handle ref_panels properly
    // For testing, create a dummy reference panel if none provided
    def ref_panels = params.ref_panels ?: [
        ['test_ref', file("NO_FILE"), file("NO_FILE")]
    ]
    
    //
    // MODULE: Minimac4 imputation
    //
    if (ref_panels && ref_panels.size() > 0) {
        // Map the phased channel to include formatted reference panel paths
        ch_phased_with_ref = ch_phased.map { meta, vcf, vcf_index ->
            // Get the first reference panel and format paths with chromosome
            def ref_panel = ref_panels[0]
            def ref_name = ref_panel[0]
            def ref_msav_template = ref_panel[1]
            def ref_vcf_template = ref_panel[2]
            
            // Replace %s with actual chromosome from meta
            def chrm = meta.contig ?: 'chr21'  // Default to chr21 for testing
            def ref_msav = file(sprintf(ref_msav_template, chrm))
            def ref_vcf = file(sprintf(ref_vcf_template, chrm))
            
            // Return the mapped data
            [meta, vcf, vcf_index, [ref_name, ref_msav, ref_vcf]]
        }
        
        IMPUTE_MINIMAC4 ( 
            ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data -> 
                [meta, vcf, vcf_index]
            },
            ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data ->
                ref_data
            }.first()  // Use the same reference for all chunks (for now)
        )
        ch_versions = ch_versions.mix(IMPUTE_MINIMAC4.out.versions)
        ch_imputed = IMPUTE_MINIMAC4.out.imputed
        ch_info = IMPUTE_MINIMAC4.out.info
    } else {
        // If no reference panels, just pass through
        ch_imputed = ch_phased
        ch_info = Channel.empty()
    }

    emit:
    imputed  = ch_imputed                  // channel: [ val(meta), path(vcf) ]
    info     = ch_info                      // channel: [ val(meta), path(info) ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}