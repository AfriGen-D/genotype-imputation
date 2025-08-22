/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PREPROCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_FILES               } from '../../modules/local/qc/check_files'
include { CHECK_CHROMOSOME          } from '../../modules/local/qc/check_chromosome'
include { CHECK_GENOME_BUILD        } from '../../modules/local/qc/check_genome_build'
include { GET_REF_CHROMOSOMES       } from '../../modules/local/qc/get_ref_chromosomes'
include { GENERATE_CHUNK_MAP        } from '../../modules/local/qc/generate_chunk_map'
include { GENERATE_REF_MAP          } from '../../modules/local/qc/generate_ref_map'
include { CHECK_OVERLAP             } from '../../modules/local/qc/check_overlap'
include { CHECK_MISMATCH            } from '../../modules/local/qc/check_mismatch'
include { MERGE_ADJACENT_CHUNKS     } from '../../modules/local/qc/merge_adjacent_chunks'
include { QC_DUPL                   } from '../../modules/local/qc/qc_dupl'
include { SPLIT_MULTI_ALLELIC       } from '../../modules/local/qc/split_multi_allelic'
include { FILTER_MIN_AC             } from '../../modules/local/qc/filter_min_ac'
// include { TARGET_QC                 } from '../../modules/local/qc/target_qc'  // Temporarily disabled - issues with QUAL scores
include { QC_SITE_MISSINGNESS       } from '../../modules/local/qc/qc_site_missingness'
include { GET_CHROMOSOME            } from '../../modules/qc'
include { check_chromosome_vcf      } from '../../modules/qc'
include { GENERATE_CHUNKS_VCF       } from '../../modules/subset_vcf'
include { SPLIT_TARGET_TO_CHUNK     } from '../../modules/subset_vcf'

workflow PREPROCESS {
    take:
    ch_input // channel: [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()
    
    //
    // MODULE: Check input files exist
    //
    CHECK_FILES ( ch_input )
    ch_versions = ch_versions.mix(CHECK_FILES.out.versions)
    
    //
    // MODULE: Check chromosome consistency
    //
    CHECK_CHROMOSOME ( CHECK_FILES.out.vcf )
    ch_versions = ch_versions.mix(CHECK_CHROMOSOME.out.versions)
    
    //
    // MODULE: Check genome build compatibility
    //
    def expected_build = params.genome_build ?: 'b38'  // Default to b38
    CHECK_GENOME_BUILD ( 
        CHECK_CHROMOSOME.out.vcf,
        expected_build
    )
    ch_versions = ch_versions.mix(CHECK_GENOME_BUILD.out.versions)
    
    // Filter out datasets with incompatible genome builds
    // We need to match the build status files with the VCF outputs
    ch_build_status = CHECK_GENOME_BUILD.out.build_status
        .map { build_file ->
            // Extract dataset ID from filename
            def dataset_id = build_file.name.replace('.build_check.txt', '')
            // Read the file and extract info
            def content = build_file.text
            def status = content.contains("PASS") ? "PASS" : "FAIL"
            
            // Extract detected build info
            def detected_build = "unknown"
            def first_chrom = "unknown"
            content.eachLine { line ->
                if (line.startsWith("First chromosome found:")) {
                    first_chrom = line.replace("First chromosome found:", "").trim()
                }
                if (line.startsWith("Detected build:")) {
                    detected_build = line.split(":")[1].split("\\(")[0].trim()
                }
            }
            
            [dataset_id, status, detected_build, first_chrom]
        }
    
    ch_vcf_build_validated = CHECK_GENOME_BUILD.out.vcf
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .join(ch_build_status, by: 0)
        .filter { dataset_id, meta, vcf, status, detected_build, first_chrom ->
            if (status == "PASS") {
                log.info "✓ Dataset ${dataset_id} passed genome build check"
                log.info "  Detected: ${detected_build} (example: ${first_chrom})"
                log.info "  Expected: ${expected_build} - Compatible ✓"
                return true
            } else {
                log.warn "================================================================"
                log.warn "DATASET EXCLUDED: ${dataset_id}"
                log.warn "----------------------------------------------------------------"
                log.warn "  Reason: Genome build incompatibility"
                log.warn "  "
                log.warn "  Dataset details:"
                log.warn "    • First chromosome: ${first_chrom}"
                log.warn "    • Detected build: ${detected_build}"
                log.warn "    • Format: ${detected_build == 'b37' ? 'chromosomes without prefix (1, 2, 3, X, Y)' : 'chromosomes with chr prefix (chr1, chr2, chr3, chrX, chrY)'}"
                log.warn "  "
                log.warn "  Expected configuration:"
                log.warn "    • Required build: ${expected_build}"
                log.warn "    • Format: ${expected_build == 'b37' ? 'chromosomes without prefix (1, 2, 3, X, Y)' : 'chromosomes with chr prefix (chr1, chr2, chr3, chrX, chrY)'}"
                log.warn "  "
                log.warn "  Action taken: Dataset will be skipped from analysis"
                log.warn "  "
                log.warn "  Solution: Provide datasets that match the ${expected_build} genome build"
                log.warn "           OR change 'genome_build' parameter in config to '${detected_build}'"
                log.warn "================================================================"
                return false
            }
        }
        .map { dataset_id, meta, vcf, status, detected_build, first_chrom -> [meta, vcf] }
    
    // Check if we have any valid datasets left
    ch_vcf_build_validated
        .count()
        .subscribe { count ->
            if (count == 0) {
                log.error "================================================================"
                log.error "ERROR: No datasets passed genome build validation!"
                log.error "Expected genome build: ${expected_build}"
                log.error "Please check that your input datasets match the expected build:"
                log.error "  - b37/hg19: chromosomes as '1', '2', '3', etc."
                log.error "  - b38/hg38: chromosomes as 'chr1', 'chr2', 'chr3', etc."
                log.error "================================================================"
                error "Pipeline stopped: No compatible datasets found"
            }
        }
    
    //
    // MODULE: Generate chunks for parallel processing FIRST
    //
    ch_vcf_for_map = ch_vcf_build_validated
        .map{ meta, vcf -> 
            [meta.id, vcf]
        }
    
    GET_CHROMOSOME(ch_vcf_for_map)
    
    ch_chr_for_chunks = GET_CHROMOSOME.out
        .map{ dataset, vcf, map_file -> 
            check_chromosome_vcf(dataset, vcf, map_file, params.chromosomes ?: '')
        }
        .map{ dataset, vcf, map_file, chrms -> 
            [dataset, vcf, map_file, chrms.unique().join(','), params.chunk_size ?: 5000000]
        }
    
    GENERATE_CHUNKS_VCF(ch_chr_for_chunks)
    
    ch_chunks = GENERATE_CHUNKS_VCF.out.flatMap{ dataset, vcf, chunk_file ->
        file(chunk_file).readLines().collect{ chunk_data ->
            def parts = chunk_data.trim().split(',')
            if (parts.size() >= 3) {
                def chrm = parts[0]
                def start = parts[1]
                def end = parts[2]
                def tagName = "${chrm}_${start}_${end}"
                [dataset, chrm, start, end, tagName, vcf]
            }
        }.findAll { it != null }
    }
    
    SPLIT_TARGET_TO_CHUNK(ch_chunks)
    
    ch_vcf_chunks = SPLIT_TARGET_TO_CHUNK.out
        .map{ dataset, chrm, start, end, tagName, vcf_chunk, vcf_index ->
            def meta = [:]
            meta.id = "${dataset}_${chrm}_${start}_${end}"  // Include chunk info in ID
            meta.sample = dataset  // Keep original sample name
            meta.chunk = "${chrm}_${start}_${end}"
            meta.contig = chrm
            meta.start = start
            meta.end = end
            [meta, vcf_chunk]
        }
    
    //
    // MODULE: Check overlap between chunks and reference panel
    //
    def ref_panels = params.ref_panels ?: []
    if (ref_panels && ref_panels.size() > 0) {
        def ref_panel = ref_panels[0]
        
        // First, get the list of chromosomes available in the reference panel
        // This will also validate that the reference panel build matches expected build
        GET_REF_CHROMOSOMES ( 
            ref_panel,
            expected_build 
        )
        ch_versions = ch_versions.mix(GET_REF_CHROMOSOMES.out.versions)
        
        // Read the chromosome list
        ch_ref_chroms = GET_REF_CHROMOSOMES.out.chromosomes
            .splitText()
            .map { it.trim() }
            .collect()
            .map { chroms -> chroms as Set }
        
        // Filter chunks to only those with matching chromosomes in reference panel
        ch_vcf_chunks_filtered_by_chrom = ch_vcf_chunks
            .combine(ch_ref_chroms)
            .map { meta, vcf, ref_chroms ->
                def chrm = meta.contig
                if (ref_chroms.contains(chrm)) {
                    [meta, vcf]
                } else {
                    // Check if this is a naming convention issue
                    def alt_chrm = chrm.startsWith('chr') ? chrm.substring(3) : "chr${chrm}"
                    if (ref_chroms.contains(alt_chrm)) {
                        log.warn "================================================================"
                        log.warn "CHUNK SKIPPED: ${meta.id}"
                        log.warn "----------------------------------------------------------------"
                        log.warn "  Reason: Chromosome naming mismatch"
                        log.warn "  "
                        log.warn "  Chunk chromosome: '${chrm}'"
                        log.warn "  Reference panel has: '${alt_chrm}' (but not '${chrm}')"
                        log.warn "  "
                        log.warn "  This indicates a genome build mismatch:"
                        log.warn "    • ${chrm.startsWith('chr') ? 'Chunk uses b38 format' : 'Chunk uses b37 format'}"
                        log.warn "    • ${alt_chrm.startsWith('chr') ? 'Reference needs b38 format' : 'Reference needs b37 format'}"
                        log.warn "  "
                        log.warn "  Action: This chunk will be excluded from processing"
                        log.warn "================================================================"
                    } else {
                        log.warn "================================================================"
                        log.warn "CHUNK SKIPPED: ${meta.id}"
                        log.warn "  Chromosome '${chrm}' not found in reference panel"
                        log.warn "  Available chromosomes: ${ref_chroms.take(5).join(', ')}..."
                        log.warn "================================================================"
                    }
                    null
                }
            }
            .filter { it != null }
        
        // Generate map for chunks that passed chromosome filter
        GENERATE_CHUNK_MAP ( ch_vcf_chunks_filtered_by_chrom )
        ch_versions = ch_versions.mix(GENERATE_CHUNK_MAP.out.versions)
        
        // Generate map for reference panel (with chunk regions)
        ch_ref_input = ch_vcf_chunks_filtered_by_chrom.map { meta, vcf ->
            def chrm = meta.contig
            def ref_bcf_path = sprintf(ref_panel[2], chrm)
            def ref_panel_bcf = file(ref_bcf_path)
            def ref_panel_index = file("${ref_bcf_path}.csi")
            [meta, ref_panel_bcf, ref_panel_index]
        }
        GENERATE_REF_MAP ( ch_ref_input )
        ch_versions = ch_versions.mix(GENERATE_REF_MAP.out.versions)
        
        // Combine maps and check overlap
        ch_maps = GENERATE_CHUNK_MAP.out.map
            .join(GENERATE_REF_MAP.out.map)
            .map { meta, chunk_map, ref_map ->
                [meta, chunk_map, ref_map]
            }
        
        CHECK_OVERLAP ( ch_maps )
        ch_versions = ch_versions.mix(CHECK_OVERLAP.out.versions)
        
        // Separate passed and failed chunks based on overlap status
        ch_overlap_results = ch_vcf_chunks_filtered_by_chrom
            .join(CHECK_OVERLAP.out.overlap)
            .branch { meta, vcf, overlap_txt, status_file ->
                def status = status_file.text.trim()
                passed: status == "PASS"
                    return [meta, vcf]
                failed: status == "FAIL"
                    return [meta, vcf]
            }
        
        // Log the failed chunks
        ch_overlap_results.failed.subscribe { meta, vcf ->
            log.warn "Chunk ${meta.id} failed overlap check - will attempt to merge with adjacent chunk"
        }
        
        // For failed chunks, find adjacent passing chunks to merge with
        ch_chunks_to_merge = ch_overlap_results.failed
            .map { meta, vcf -> 
                // Create a key based on sample and chromosome
                def key = "${meta.sample}_${meta.contig}"
                [key, meta, vcf]
            }
            .combine(
                ch_overlap_results.passed.map { meta, vcf ->
                    def key = "${meta.sample}_${meta.contig}"
                    [key, meta, vcf]
                },
                by: 0  // Combine by sample_chromosome key
            )
            .map { key, failed_meta, failed_vcf, passed_meta, passed_vcf ->
                // Calculate distance between chunks
                def distance = Math.abs(failed_meta.start.toLong() - passed_meta.start.toLong())
                [passed_meta, failed_meta, passed_vcf, failed_vcf, distance]
            }
            .groupTuple(by: [0, 1])  // Group by passed and failed meta
            .map { passed_meta, failed_meta, passed_vcf, failed_vcf, distances ->
                // Find the minimum distance (closest chunk)
                def min_idx = distances.indexOf(distances.min())
                [passed_meta, failed_meta, passed_vcf[min_idx], failed_vcf[min_idx]]
            }
        
        // Merge failed chunks with their nearest passing chunks
        MERGE_ADJACENT_CHUNKS ( ch_chunks_to_merge )
        ch_versions = ch_versions.mix(MERGE_ADJACENT_CHUNKS.out.versions)
        
        // Combine the merged chunks with the originally passed chunks that didn't need merging
        // Simple solution: If no merges needed, use passed chunks; otherwise filter and mix
        ch_vcf_chunks_filtered = MERGE_ADJACENT_CHUNKS.out.vcf
            .map { meta, vcf, tbi -> [meta, vcf] }
            .mix(ch_overlap_results.passed)
            .unique { meta, vcf -> meta.id }  // Keep unique chunks by ID
            .ifEmpty(ch_overlap_results.passed)  // If no merges, use all passed chunks
        
        // Check allele mismatch for chunks that passed overlap check
        ch_mismatch_input = ch_vcf_chunks_filtered.map { meta, vcf ->
            def chrm = meta.contig
            def ref_bcf_path = sprintf(ref_panel[2], chrm)
            def ref_panel_bcf = file(ref_bcf_path)
            def ref_panel_index = file("${ref_bcf_path}.csi")
            [meta, vcf, ref_panel_bcf, ref_panel_index]
        }
        
        CHECK_MISMATCH ( ch_mismatch_input )
        ch_versions = ch_versions.mix(CHECK_MISMATCH.out.versions)
        
        // Filter chunks based on mismatch check
        ch_vcf_chunks_final = ch_vcf_chunks_filtered
            .join(CHECK_MISMATCH.out.mismatch)
            .filter { meta, vcf, mismatch_txt, mismatch_status ->
                def status = mismatch_status.text.trim()
                if (status == "PASS") {
                    return true
                } else {
                    log.warn "Skipping chunk ${meta.id} due to high allele mismatch with reference panel"
                    return false
                }
            }
            .map { meta, vcf, mismatch_txt, mismatch_status ->
                [meta, vcf]
            }
        
        ch_vcf_chunks = ch_vcf_chunks_final
    } else {
        // If no reference panels configured, use all chunks as-is
        ch_vcf_chunks = ch_vcf_chunks
    }
    
    //
    // MODULE: QC for duplicates - on chunks
    //
    QC_DUPL ( ch_vcf_chunks )
    ch_versions = ch_versions.mix(QC_DUPL.out.versions)
    
    //
    // MODULE: Split multi-allelic variants
    //
    SPLIT_MULTI_ALLELIC ( QC_DUPL.out.vcf )
    ch_versions = ch_versions.mix(SPLIT_MULTI_ALLELIC.out.versions)
    
    //
    // MODULE: Filter by minimum allele count
    //
    FILTER_MIN_AC ( SPLIT_MULTI_ALLELIC.out.vcf )
    ch_versions = ch_versions.mix(FILTER_MIN_AC.out.versions)
    
    //
    // MODULE: Target QC - SKIPPED (commented out for debugging)
    // TARGET_QC module was causing issues with array genotyping data that lacks QUAL scores
    // Skipping this step to allow pipeline to proceed
    //
    // TARGET_QC ( FILTER_MIN_AC.out.vcf )
    // ch_versions = ch_versions.mix(TARGET_QC.out.versions)
    
    //
    // MODULE: Site missingness QC
    // Now takes input directly from FILTER_MIN_AC instead of TARGET_QC
    //
    QC_SITE_MISSINGNESS ( FILTER_MIN_AC.out.vcf )
    ch_versions = ch_versions.mix(QC_SITE_MISSINGNESS.out.versions)
    
    // SITES_ONLY process removed - not currently used by downstream processes

    emit:
    vcf      = QC_SITE_MISSINGNESS.out.vcf  // channel: [ val(meta), path(vcf) ] - QC'd chunks with genotypes
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}