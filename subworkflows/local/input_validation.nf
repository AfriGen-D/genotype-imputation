/*
========================================================================================
    SUBWORKFLOW: INPUT_VALIDATION
========================================================================================
*/

include { CHECK_VCF_FORMAT } from '../../modules/local/qc/check_vcf_format'
include { CHECK_REFERENCE_COMPATIBILITY } from '../../modules/local/qc/check_reference'
include { VALIDATE_CHROMOSOMES } from '../../modules/local/qc/validate_chromosomes'
include { CHECK_SAMPLE_OVERLAP } from '../../modules/local/qc/check_samples'

workflow INPUT_VALIDATION {
    take:
    ch_input_vcfs      // channel: [ val(meta), path(vcf) ]
    ch_reference_panels // channel: [ val(ref_name), path(vcf), path(m3vcf) ]
    genome_build       // val: 'b37' or 'b38'

    main:
    ch_versions = Channel.empty()
    
    // Validate VCF format and integrity
    CHECK_VCF_FORMAT(ch_input_vcfs)
    ch_versions = ch_versions.mix(CHECK_VCF_FORMAT.out.versions)
    
    // Check chromosome naming consistency
    VALIDATE_CHROMOSOMES(
        CHECK_VCF_FORMAT.out.vcf,
        genome_build
    )
    ch_versions = ch_versions.mix(VALIDATE_CHROMOSOMES.out.versions)
    
    // Verify reference panel compatibility
    CHECK_REFERENCE_COMPATIBILITY(
        CHECK_VCF_FORMAT.out.vcf,
        ch_reference_panels
    )
    ch_versions = ch_versions.mix(CHECK_REFERENCE_COMPATIBILITY.out.versions)
    
    // Check sample overlap between input and reference
    CHECK_SAMPLE_OVERLAP(
        CHECK_VCF_FORMAT.out.vcf,
        ch_reference_panels
    )
    ch_versions = ch_versions.mix(CHECK_SAMPLE_OVERLAP.out.versions)
    
    emit:
    vcf      = CHECK_VCF_FORMAT.out.vcf
    stats    = CHECK_SAMPLE_OVERLAP.out.report
    versions = ch_versions
}