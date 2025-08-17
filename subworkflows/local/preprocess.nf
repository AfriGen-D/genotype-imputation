/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PREPROCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_FILES               } from '../../modules/local/qc/check_files'
include { CHECK_CHROMOSOME          } from '../../modules/local/qc/check_chromosome'
include { GENERATE_CHUNK_MAP        } from '../../modules/local/qc/generate_chunk_map'
include { GENERATE_REF_MAP          } from '../../modules/local/qc/generate_ref_map'
include { CHECK_OVERLAP             } from '../../modules/local/qc/check_overlap'
include { CHECK_MISMATCH            } from '../../modules/local/qc/check_mismatch'
include { QC_DUPL                   } from '../../modules/local/qc/qc_dupl'
include { SPLIT_MULTI_ALLELIC       } from '../../modules/local/qc/split_multi_allelic'
include { FILTER_MIN_AC             } from '../../modules/local/qc/filter_min_ac'
include { TARGET_QC                 } from '../../modules/local/qc/target_qc'
include { QC_SITE_MISSINGNESS       } from '../../modules/local/qc/qc_site_missingness'
include { SITES_ONLY                } from '../../modules/local/qc/sites_only'
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
    // MODULE: Generate chunks for parallel processing FIRST
    //
    ch_vcf_for_map = CHECK_CHROMOSOME.out.vcf
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
        .map{ dataset, chrm, start, end, tagName, vcf_chunk ->
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
        def ref_bcf_path = sprintf(ref_panel[2], 'chr21')
        ref_panel_bcf = file(ref_bcf_path)
        ref_panel_index = file("${ref_bcf_path}.csi")
        
        // Generate map for chunks
        GENERATE_CHUNK_MAP ( ch_vcf_chunks )
        ch_versions = ch_versions.mix(GENERATE_CHUNK_MAP.out.versions)
        
        // Generate map for reference panel (with chunk regions)
        ch_ref_input = ch_vcf_chunks.map { meta, vcf ->
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
        
        // Filter chunks based on overlap status - only pass chunks that meet threshold
        ch_vcf_chunks_filtered = ch_vcf_chunks
            .join(CHECK_OVERLAP.out.overlap)
            .filter { meta, vcf, overlap_txt, status_file ->
                def status = status_file.text.trim()
                if (status == "PASS") {
                    return true
                } else {
                    log.warn "Skipping chunk ${meta.id} due to insufficient overlap with reference panel"
                    return false
                }
            }
            .map { meta, vcf, overlap_txt, status_file ->
                [meta, vcf]
            }
        
        // Check allele mismatch for chunks that passed overlap check
        ch_mismatch_input = ch_vcf_chunks_filtered.map { meta, vcf ->
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
    // MODULE: Target QC
    //
    TARGET_QC ( FILTER_MIN_AC.out.vcf )
    ch_versions = ch_versions.mix(TARGET_QC.out.versions)
    
    //
    // MODULE: Site missingness QC
    //
    QC_SITE_MISSINGNESS ( TARGET_QC.out.vcf )
    ch_versions = ch_versions.mix(QC_SITE_MISSINGNESS.out.versions)
    
    //
    // MODULE: Extract sites only (for reference, but not for phasing)
    //
    SITES_ONLY ( QC_SITE_MISSINGNESS.out.vcf )
    ch_versions = ch_versions.mix(SITES_ONLY.out.versions)

    emit:
    vcf      = QC_SITE_MISSINGNESS.out.vcf  // channel: [ val(meta), path(vcf) ] - QC'd chunks with genotypes
    sites    = SITES_ONLY.out.sites         // channel: [ val(meta), path(sites) ]
    sites_only = SITES_ONLY.out.vcf         // channel: [ val(meta), path(vcf) ] - Sites only (no genotypes)
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}