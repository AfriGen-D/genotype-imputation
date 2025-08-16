/*
========================================================================================
    SUBWORKFLOW: PHASING
========================================================================================
*/

include { CHUNK_GENOME } from '../../modules/local/vcf/chunk_genome'
include { EAGLE_PHASE } from '../../modules/local/phasing/eagle'
include { SHAPEIT_PHASE } from '../../modules/local/phasing/shapeit'
include { BEAGLE_PHASE } from '../../modules/local/phasing/beagle'
include { MERGE_PHASED_CHUNKS } from '../../modules/local/vcf/merge_chunks'

workflow PHASING {
    take:
    ch_vcfs           // channel: [ val(meta), path(vcf) ]
    ch_reference      // channel: [ val(ref_name), path(vcf) ]
    ch_genetic_map    // path: genetic map file
    phasing_tool      // val: 'eagle', 'shapeit', or 'beagle'
    chunk_size        // val: chunk size in base pairs

    main:
    ch_versions = Channel.empty()
    
    // Check if input is already chunked (has chunk_id in meta)
    // If already chunked, use directly; otherwise create chunks
    ch_vcfs
        .branch {
            chunked: it[0].containsKey('chunk_id')
            not_chunked: !it[0].containsKey('chunk_id')
        }
        .set { ch_vcf_branched }
    
    // Create chunks only for non-chunked VCFs
    CHUNK_GENOME(
        ch_vcf_branched.not_chunked,
        chunk_size
    )
    ch_versions = ch_versions.mix(CHUNK_GENOME.out.versions)
    
    // Combine already chunked and newly chunked VCFs
    ch_chunks_combined = ch_vcf_branched.chunked
        .mix(CHUNK_GENOME.out.chunks)
    
    // Combine chunks with reference panels
    ch_chunks_with_ref = ch_chunks_combined
        .combine(ch_reference)
    
    // Phase using selected tool
    if (phasing_tool == 'eagle') {
        EAGLE_PHASE(
            ch_chunks_with_ref,
            ch_genetic_map
        )
        ch_phased = EAGLE_PHASE.out.phased_vcf
        ch_versions = ch_versions.mix(EAGLE_PHASE.out.versions)
    } else if (phasing_tool == 'shapeit') {
        SHAPEIT_PHASE(
            ch_chunks_with_ref,
            ch_genetic_map
        )
        ch_phased = SHAPEIT_PHASE.out.phased_vcf
        ch_versions = ch_versions.mix(SHAPEIT_PHASE.out.versions)
    } else if (phasing_tool == 'beagle') {
        BEAGLE_PHASE(
            ch_chunks_with_ref,
            ch_genetic_map
        )
        ch_phased = BEAGLE_PHASE.out.phased_vcf
        ch_versions = ch_versions.mix(BEAGLE_PHASE.out.versions)
    } else {
        error "Unknown phasing tool: ${phasing_tool}"
    }
    
    // Group phased chunks by sample
    ch_grouped_chunks = ch_phased
        .map { meta, vcf, ref_name, ref_vcf -> 
            [ meta.id, meta, vcf ]
        }
        .groupTuple(by: 0)
        .map { id, metas, vcfs -> 
            [ metas[0], vcfs ]
        }
    
    // Merge phased chunks back together
    MERGE_PHASED_CHUNKS(ch_grouped_chunks)
    ch_versions = ch_versions.mix(MERGE_PHASED_CHUNKS.out.versions)
    
    emit:
    phased_vcf = MERGE_PHASED_CHUNKS.out.merged_vcf
    versions   = ch_versions
}