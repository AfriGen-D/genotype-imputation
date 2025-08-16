/*
========================================================================================
    SUBWORKFLOW: VCF_CHUNKING
    Split VCF files into chunks for parallel processing
========================================================================================
*/

include { CREATE_CHUNK_LIST } from '../../modules/local/vcf/create_chunk_list'
include { EXTRACT_VCF_CHUNK } from '../../modules/local/vcf/extract_vcf_chunk'
include { CHECK_CHUNK_OVERLAP } from '../../modules/local/qc/check_chunk_overlap'
include { AGGREGATE_OVERLAP_REPORTS } from '../../modules/local/qc/aggregate_overlap_reports'
include { MERGE_LOW_OVERLAP_CHUNKS } from '../../modules/local/vcf/merge_low_overlap_chunks'

workflow VCF_CHUNKING {
    take:
    ch_vcfs              // channel: [ val(meta), path(vcf), path(vcf_index) ]
    ch_reference_panels  // channel: [ val(ref_name), path(m3vcf), path(vcf) ]
    chunk_size           // val: genomic region size in base pairs

    main:
    ch_versions = Channel.empty()
    
    // Step 1: Create chunk list (fast, single process)
    CREATE_CHUNK_LIST(
        ch_vcfs,
        chunk_size
    )
    ch_versions = ch_versions.mix(CREATE_CHUNK_LIST.out.versions)
    
    // Step 2: Parse chunk list and create channel for parallel extraction
    ch_chunks = CREATE_CHUNK_LIST.out.chunk_list
        .map { meta, chunk_file ->
            // Read the chunk file and create entries for each chunk
            def chunks = []
            chunk_file.eachLine { line ->
                if (!line.startsWith('chunk_name')) {  // Skip header
                    def parts = line.split('\t')
                    if (parts.size() >= 2) {
                        chunks.add([meta, parts[0], parts[1]])  // meta, chunk_name, region
                    }
                }
            }
            return chunks
        }
        .flatMap()  // Flatten the list of chunks
        .combine(ch_vcfs, by: 0)  // Combine with original VCF by meta
        .map { meta, chunk_name, region, vcf, vcf_index ->
            [ meta, chunk_name, region, vcf, vcf_index ]
        }
    
    // Step 3: Extract chunks in parallel
    EXTRACT_VCF_CHUNK(ch_chunks)
    ch_versions = ch_versions.mix(EXTRACT_VCF_CHUNK.out.versions)
    
    // Step 4: Check overlap between target chunks and reference panels
    // Combine each extracted chunk with reference panels for overlap checking
    ch_overlap_input = ch_chunks
        .map { meta, chunk_name, region, vcf, vcf_index ->
            [ meta, chunk_name, region ]
        }
        .join(
            EXTRACT_VCF_CHUNK.out.chunk
                .map { meta, vcf, index ->
                    def chunk_name = vcf.baseName
                    [ meta, chunk_name, vcf, index ]
                },
            by: [0, 1]  // Join by meta and chunk_name
        )
        .combine(ch_reference_panels)  // Combine with all reference panels
        .map { meta, chunk_name, region, target_vcf, target_index, ref_name, m3vcf, ref_vcf ->
            [ meta, chunk_name, region, target_vcf, target_index, ref_name, ref_vcf ]
        }
    
    CHECK_CHUNK_OVERLAP(ch_overlap_input)
    ch_versions = ch_versions.mix(CHECK_CHUNK_OVERLAP.out.versions)
    
    // Step 5: Aggregate all overlap reports
    AGGREGATE_OVERLAP_REPORTS(
        CHECK_CHUNK_OVERLAP.out.overlap_report.map { meta, chunk, ref, report -> report }.collect(),
        CHECK_CHUNK_OVERLAP.out.summary.map { meta, chunk, ref, summary -> summary }.collect()
    )
    ch_versions = ch_versions.mix(AGGREGATE_OVERLAP_REPORTS.out.versions)
    
    // Step 6: Merge chunks with low overlap
    // Group chunks by sample with their overlap matrix
    ch_chunks_for_merge = EXTRACT_VCF_CHUNK.out.chunk
        .map { meta, vcf, index -> 
            [ meta.id, meta, vcf, index ]
        }
        .groupTuple(by: 0)
        .map { id, metas, vcfs, indices -> 
            [ metas[0], vcfs, indices ]
        }
        .join(
            AGGREGATE_OVERLAP_REPORTS.out.matrix
                .map { matrix -> 
                    // Extract sample ID from matrix filename if needed
                    [ "default", matrix ]  // Using default key for now
                },
            by: 0
        )
        .map { id, meta, vcfs, indices, matrix ->
            [ meta, vcfs, indices, matrix ]
        }
    
    // Merge chunks with overlap below threshold
    MERGE_LOW_OVERLAP_CHUNKS(
        ch_chunks_for_merge,
        params.overlap_min_ratio ?: 0.00001  // Use overlap_min_ratio parameter (0.001% default)
    )
    ch_versions = ch_versions.mix(MERGE_LOW_OVERLAP_CHUNKS.out.versions)
    
    // Use merged chunks for downstream processing
    ch_final_chunks = MERGE_LOW_OVERLAP_CHUNKS.out.merged_chunks
        .transpose()
        .map { meta, vcf, index ->
            // Add chunk info to meta
            def chunk_id = vcf.baseName.replaceAll(/.*_merged_chunk_/, '')
            def new_meta = meta + [chunk_id: chunk_id, merged: true]
            [ new_meta, vcf, index ]
        }
    
    emit:
    chunks            = ch_final_chunks
    overlap_summary   = AGGREGATE_OVERLAP_REPORTS.out.summary
    overlap_matrix    = AGGREGATE_OVERLAP_REPORTS.out.matrix
    merge_report      = MERGE_LOW_OVERLAP_CHUNKS.out.merge_report
    stats            = EXTRACT_VCF_CHUNK.out.stats.collect()
    summary          = CREATE_CHUNK_LIST.out.summary
    versions         = ch_versions
}