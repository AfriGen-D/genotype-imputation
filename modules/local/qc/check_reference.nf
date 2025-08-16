process CHECK_REFERENCE_COMPATIBILITY {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf)
    tuple val(ref_name), path(ref_vcf), path(ref_m3vcf)
    
    output:
    tuple val(meta), path("${meta.id}_ref_check.txt"), emit: report
    tuple val(meta), path("${meta.id}_ref_timing.log"), emit: timing
    path "versions.yml", emit: versions
    
    script:
    """
    # Initialize timing log
    echo "=== CHECK_REFERENCE_COMPATIBILITY Timing Report ===" > ${meta.id}_ref_timing.log
    echo "Process started at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${meta.id}_ref_timing.log
    START_TIME=\$(date +%s)
    
    # Check chromosome naming convention
    echo "Step 1: Checking chromosome formats started at \$(date '+%H:%M:%S')" >> ${meta.id}_ref_timing.log
    STEP_START=\$(date +%s)
    bcftools query -f '%CHROM\\n' $vcf | head -1 > chr_format.txt
    bcftools query -f '%CHROM\\n' $ref_vcf | head -1 > ref_chr_format.txt
    echo "  Chromosome format check completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_ref_timing.log
    
    # Generate compatibility report
    echo "Step 2: Generating compatibility report started at \$(date '+%H:%M:%S')" >> ${meta.id}_ref_timing.log
    STEP_START=\$(date +%s)
    echo "Reference Compatibility Check for ${meta.id}" > ${meta.id}_ref_check.txt
    echo "VCF chromosome format: \$(cat chr_format.txt)" >> ${meta.id}_ref_check.txt
    echo "Reference chromosome format: \$(cat ref_chr_format.txt)" >> ${meta.id}_ref_check.txt
    
    # Check if formats match
    if grep -q "chr" chr_format.txt && ! grep -q "chr" ref_chr_format.txt; then
        echo "WARNING: Chromosome naming mismatch detected" >> ${meta.id}_ref_check.txt
    elif ! grep -q "chr" chr_format.txt && grep -q "chr" ref_chr_format.txt; then
        echo "WARNING: Chromosome naming mismatch detected" >> ${meta.id}_ref_check.txt
    else
        echo "Chromosome naming is compatible" >> ${meta.id}_ref_check.txt
    fi
    echo "  Compatibility report completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_ref_timing.log
    
    # Final timing report
    END_TIME=\$(date +%s)
    echo "Process completed at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${meta.id}_ref_timing.log
    echo "Total execution time: \$((END_TIME - START_TIME)) seconds" >> ${meta.id}_ref_timing.log
    echo "=== End of Timing Report ===" >> ${meta.id}_ref_timing.log
    
    # Display timing summary
    cat ${meta.id}_ref_timing.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}