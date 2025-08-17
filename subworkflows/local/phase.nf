/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PHASE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENERATE_CHUNKS_VCF       } from '../../modules/local/subset_vcf/generate_chunks_vcf'
// include { SPLIT_TARGET_TO_CHUNK     } from '../../modules/local/subset_vcf/split_target_to_chunk'
include { EAGLE_PHASING             } from '../../modules/local/phasing/eagle_phasing'

workflow PHASE {
    take:
    ch_vcf // channel: [ val(meta), path(vcf), path(vcf_index) ]

    main:
    ch_versions = Channel.empty()
    
    // Phasing step - ch_vcf already has index from QC_SITE_MISSINGNESS output (with genotypes)
    ch_vcf_with_index = ch_vcf
    
    //
    // MODULE: Eagle phasing with reference panel
    //
    // Get genetic map
    def genetic_map = params.eagle_genetic_map ? file(params.eagle_genetic_map) : []
    
    // Get reference panel BCF for chr21 (since we're testing with chr21 data)
    // TODO: Make this dynamic based on actual chromosome being processed
    def ref_panel_bcf = []
    def ref_panel_index = []
    if (params.ref_panels && params.ref_panels.size() > 0) {
        // Get the first reference panel and format for chr21
        def ref_panel = params.ref_panels[0]
        def ref_bcf_path = sprintf(ref_panel[2], 'chr21')  // Use chr21 for now
        ref_panel_bcf = file(ref_bcf_path)
        ref_panel_index = file("${ref_bcf_path}.csi")
    }
    
    EAGLE_PHASING ( 
        ch_vcf_with_index,
        genetic_map,
        ref_panel_bcf,
        ref_panel_index
    )
    ch_versions = ch_versions.mix(EAGLE_PHASING.out.versions)

    emit:
    phased   = EAGLE_PHASING.out.phased  // channel: [ val(meta), path(vcf) ]
    versions = ch_versions                // channel: [ path(versions.yml) ]
}