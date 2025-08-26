/*
========================================================================================
    SUBWORKFLOW: INPUT_VALIDATION
========================================================================================
*/

include { CHECK_VCF_FORMAT } from '../../modules/local/qc/check_vcf_format'
include { CHECK_REFERENCE_COMPATIBILITY } from '../../modules/local/qc/check_reference'
include { VALIDATE_CHROMOSOMES } from '../../modules/local/qc/validate_chromosomes'
include { CHECK_SAMPLE_OVERLAP } from '../../modules/local/qc/check_samples'
include { VCF_CHUNKING } from '../local/vcf_chunking'

workflow INPUT_VALIDATION {
    take:
    ch_input_vcfs      // channel: [ val(meta), path(vcf) ]
    ch_reference_panels // channel: [ val(ref_name), path(vcf), path(m3vcf) ]
    genome_build       // val: 'b37' or 'b38'

    main:
    ch_versions = Channel.empty()
    ch_timing_logs = Channel.empty()
    
    // Always validate VCF format and integrity (critical)
    CHECK_VCF_FORMAT(ch_input_vcfs)
    ch_versions = ch_versions.mix(CHECK_VCF_FORMAT.out.versions)
    ch_timing_logs = ch_timing_logs.mix(CHECK_VCF_FORMAT.out.timing)
    
    // Optionally chunk VCF for parallel processing
    if (params.enable_chunking && params.chunk_size > 0) {
        // Prepare VCF with index for chunking
        ch_vcf_with_index = CHECK_VCF_FORMAT.out.vcf
            .map { meta, vcf -> 
                def vcf_index = file("${vcf}.tbi")
                [meta, vcf, vcf_index]
            }
        
        // Perform chunking with overlap checking
        VCF_CHUNKING(
            ch_vcf_with_index,
            ch_reference_panels,
            params.chunk_size
        )
        ch_versions = ch_versions.mix(VCF_CHUNKING.out.versions)
        
        // Use chunks for output
        ch_vcf_output = VCF_CHUNKING.out.chunks
            .map { meta, vcf, index ->
                [meta, vcf]
            }
    } else {
        // Process the whole VCF directly (no chunking)
        ch_vcf_output = CHECK_VCF_FORMAT.out.vcf
    }
    
    // Optional: Check chromosome naming consistency
    if (!params.skip_chr_validation) {
        VALIDATE_CHROMOSOMES(
            ch_vcf_output,
            genome_build
        )
        ch_versions = ch_versions.mix(VALIDATE_CHROMOSOMES.out.versions)
        ch_timing_logs = ch_timing_logs.mix(VALIDATE_CHROMOSOMES.out.timing)
    }
    
    // Note: Reference compatibility and sample overlap checks are commented out
    // in the original workflow, so keeping them commented but with timing support
    
    // Collect and display all timing logs
    ch_timing_logs
        .collectFile(name: 'validation_timing_summary.log', storeDir: "${params.outdir}/logs")
        .view { "Validation timing logs saved to: ${params.outdir}/logs/validation_timing_summary.log" }
    
    emit:
    vcf      = ch_vcf_output  // Emit validated VCF
    timing   = ch_timing_logs
    versions = ch_versions
}