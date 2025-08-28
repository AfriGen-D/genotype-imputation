/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: IMPUTE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IMPUTE_MINIMAC4           } from '../../modules/local/impute/impute_minimac4'

workflow IMPUTE {
    take:
    ch_phased // channel: [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()
    ch_ref_vcf = Channel.empty()
    
    // Simplified imputation - handle ref_panels properly
    // For testing, create a dummy reference panel if none provided
    def ref_panels = params.ref_panels ?: [
        ['test_ref', file("NO_FILE"), file("NO_FILE")]
    ]
    
    //
    // MODULE: Minimac4 imputation - Standard version only
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
            def chrm = meta.contig
            def ref_msav = file(ref_msav_template.replace('%s', chrm))
            def ref_vcf = file(ref_vcf_template.replace('%s', chrm))
            
            // Return the mapped data
            [meta, vcf, vcf_index, [ref_name, ref_msav, ref_vcf]]
        }
        
        // Run standard imputation
        IMPUTE_MINIMAC4 ( 
            ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data -> 
                [meta, vcf, vcf_index]
            },
            ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data ->
                [ref_data[0], ref_data[1], ref_data[2]]
            }
        )
        ch_versions = ch_versions.mix(IMPUTE_MINIMAC4.out.versions)
        ch_imputed = IMPUTE_MINIMAC4.out.imputed
        ch_info = IMPUTE_MINIMAC4.out.info
        
        // Extract reference VCF paths for frequency comparison
        ch_ref_vcf = ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data ->
            [meta, ref_data[0], ref_data[2]]   // [ val(meta), val(ref_name), path(ref_vcf) ]
        }
    } else {
        // If no reference panels, just pass through
        ch_imputed = ch_phased
        ch_info = Channel.empty()
        ch_ref_vcf = Channel.empty()
    }

    emit:
    imputed  = ch_imputed                  // channel: [ val(meta), val(ref_name), path(vcf), path(vcf_index) ]
    info     = ch_info                      // channel: [ val(meta), path(info) ]
    ref_vcf  = ch_ref_vcf                   // channel: [ val(meta), val(ref_name), path(ref_vcf) ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}