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
    
    // Get reference panel BCF dynamically based on the chromosome being processed
    // Since Eagle processes each chunk separately, we need to pass the correct reference panel for each
    // Create an augmented channel that includes reference panel paths
    ch_eagle_input = ch_vcf_with_index.map { meta, vcf, index ->
        if (params.ref_panels && params.ref_panels.size() > 0) {
            def ref_panel = params.ref_panels[0]
            def chrm = meta.contig
            def ref_bcf_path = sprintf(ref_panel[2], chrm)
            def ref_panel_bcf = file(ref_bcf_path)
            def ref_panel_index = file("${ref_bcf_path}.csi")
            [meta, vcf, index, ref_panel_bcf, ref_panel_index]
        } else {
            [meta, vcf, index, file("NO_FILE"), file("NO_FILE")]
        }
    }
    
    // Since EAGLE_PHASING expects separate inputs, we need to restructure
    // Extract each component for the module input
    // NOTE: We must NOT use .first() on reference panels as each chromosome needs its own reference
    EAGLE_PHASING (
        ch_eagle_input.map { meta, vcf, index, ref_bcf, ref_idx -> [meta, vcf, index] },
        genetic_map,
        ch_eagle_input.map { meta, vcf, index, ref_bcf, ref_idx -> ref_bcf },
        ch_eagle_input.map { meta, vcf, index, ref_bcf, ref_idx -> ref_idx }
    )
    ch_versions = ch_versions.mix(EAGLE_PHASING.out.versions)

    emit:
    phased   = EAGLE_PHASING.out.phased  // channel: [ val(meta), path(vcf) ]
    versions = ch_versions                // channel: [ path(versions.yml) ]
}