/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: IMPUTE_SAFE
    Safe imputation with pre-validation to handle empty chunks gracefully
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VALIDATE_IMPUTATION_CHUNK } from '../../modules/local/impute/validate_imputation_chunk'
include { IMPUTE_MINIMAC4_SAFE      } from '../../modules/local/impute/impute_minimac4_safe'

workflow IMPUTE_SAFE {
    take:
    ch_phased // channel: [ val(meta), path(vcf), path(vcf_index) ]

    main:
    ch_versions = Channel.empty()
    ch_ref_vcf = Channel.empty()
    ch_skipped_chunks = Channel.empty()
    
    // Get reference panels
    def ref_panels = params.ref_panels ?: [
        ['test_ref', file("NO_FILE"), file("NO_FILE")]
    ]
    
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
            def ref_msav = file(sprintf(ref_msav_template, chrm))
            def ref_vcf = file(sprintf(ref_vcf_template, chrm))
            
            // Return the mapped data
            [meta, vcf, vcf_index, [ref_name, ref_msav, ref_vcf]]
        }
        
        // Optional: First validate chunks if enabled
        if (params.validate_chunks == true) {
            VALIDATE_IMPUTATION_CHUNK(
                ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data -> 
                    [meta, vcf, vcf_index]
                },
                ch_phased_with_ref.map { meta, vcf, vcf_index, ref_data ->
                    [ref_data[0], ref_data[1], ref_data[2]]
                }
            )
            ch_versions = ch_versions.mix(VALIDATE_IMPUTATION_CHUNK.out.versions)
            
            // Filter chunks based on validation status
            ch_validated = VALIDATE_IMPUTATION_CHUNK.out.status
                .join(ch_phased_with_ref)
                .filter { meta, status_file, vcf, vcf_index, ref_data ->
                    def status = status_file.text.trim()
                    if (status == "FAIL") {
                        log.warn "Skipping chunk ${meta.id} due to validation failure"
                        false
                    } else {
                        true
                    }
                }
                .map { meta, status_file, vcf, vcf_index, ref_data ->
                    [meta, vcf, vcf_index, ref_data]
                }
        } else {
            // Skip validation, use all chunks
            ch_validated = ch_phased_with_ref
        }
        
        // Run safe imputation
        IMPUTE_MINIMAC4_SAFE(
            ch_validated.map { meta, vcf, vcf_index, ref_data -> 
                [meta, vcf, vcf_index]
            },
            ch_validated.map { meta, vcf, vcf_index, ref_data ->
                [ref_data[0], ref_data[1], ref_data[2]]
            }
        )
        
        ch_versions = ch_versions.mix(IMPUTE_MINIMAC4_SAFE.out.versions)
        
        // Separate successful imputations from skipped chunks
        ch_imputed = IMPUTE_MINIMAC4_SAFE.out.imputed
            .filter { it[2] != null }  // Filter out empty optional outputs
        
        ch_info = IMPUTE_MINIMAC4_SAFE.out.info
            .filter { it[1] != null }  // Filter out empty optional outputs
        
        ch_skipped_chunks = IMPUTE_MINIMAC4_SAFE.out.skipped
            .filter { it[1] != null }  // Only keep actual skipped chunks
        
        // Extract reference VCF paths for frequency comparison
        ch_ref_vcf = ch_validated.map { meta, vcf, vcf_index, ref_data ->
            [meta, ref_data[0], ref_data[2]]   // [ val(meta), val(ref_name), path(ref_vcf) ]
        }
        
    } else {
        // If no reference panels, just pass through
        ch_imputed = ch_phased.map { meta, vcf, vcf_index ->
            [meta, "no_ref", vcf, vcf_index]
        }
        ch_info = Channel.empty()
        ch_ref_vcf = Channel.empty()
    }
    
    // Log summary of skipped chunks
    ch_skipped_chunks
        .collect()
        .subscribe { skip_files ->
            if (skip_files.size() > 0) {
                log.info "======================================"
                log.info "Summary: ${skip_files.size()} chunks skipped due to insufficient variants"
                log.info "Skipped chunks saved to output directory"
                log.info "======================================"
            }
        }

    emit:
    imputed       = ch_imputed         // channel: [ val(meta), val(ref_name), path(vcf), path(vcf_index) ]
    info          = ch_info            // channel: [ val(meta), path(info) ]
    ref_vcf       = ch_ref_vcf         // channel: [ val(meta), val(ref_name), path(ref_vcf) ]
    skipped       = ch_skipped_chunks  // channel: [ val(meta), path(skipped_file) ]
    versions      = ch_versions        // channel: [ path(versions.yml) ]
}