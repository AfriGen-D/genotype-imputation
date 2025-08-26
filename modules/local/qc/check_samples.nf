process CHECK_SAMPLE_OVERLAP {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf)
    tuple val(ref_name), path(ref_vcf), path(ref_m3vcf)
    
    output:
    tuple val(meta), path("${meta.id}_sample_check.txt"), emit: report
    tuple val(meta), path("${meta.id}_sample_timing.log"), emit: timing
    path "versions.yml", emit: versions
    
    script:
    """
    # Initialize timing log
    echo "=== CHECK_SAMPLE_OVERLAP Timing Report ===" > ${meta.id}_sample_timing.log
    echo "Process started at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${meta.id}_sample_timing.log
    START_TIME=\$(date +%s)
    
    # Get sample lists
    echo "Step 1: Extracting sample lists started at \$(date '+%H:%M:%S')" >> ${meta.id}_sample_timing.log
    STEP_START=\$(date +%s)
    bcftools query -l $vcf > vcf_samples.txt
    bcftools query -l $ref_vcf > ref_samples.txt 2>/dev/null || touch ref_samples.txt
    echo "  Sample extraction completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_sample_timing.log
    
    # Count samples
    echo "Step 2: Counting samples started at \$(date '+%H:%M:%S')" >> ${meta.id}_sample_timing.log
    STEP_START=\$(date +%s)
    vcf_count=\$(wc -l < vcf_samples.txt)
    ref_count=\$(wc -l < ref_samples.txt)
    echo "  Sample counting completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_sample_timing.log
    
    # Check for overlapping samples
    echo "Step 3: Checking sample overlap started at \$(date '+%H:%M:%S')" >> ${meta.id}_sample_timing.log
    STEP_START=\$(date +%s)
    comm -12 <(sort vcf_samples.txt) <(sort ref_samples.txt) > overlapping_samples.txt || true
    overlap_count=\$(wc -l < overlapping_samples.txt)
    echo "  Sample overlap check completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_sample_timing.log
    
    # Generate report
    echo "Step 4: Generating report started at \$(date '+%H:%M:%S')" >> ${meta.id}_sample_timing.log
    STEP_START=\$(date +%s)
    echo "Sample Check for ${meta.id}" > ${meta.id}_sample_check.txt
    echo "VCF samples: \$vcf_count" >> ${meta.id}_sample_check.txt
    echo "Reference samples: \$ref_count" >> ${meta.id}_sample_check.txt
    echo "Overlapping samples: \$overlap_count" >> ${meta.id}_sample_check.txt
    
    if [ \$overlap_count -gt 0 ]; then
        echo "WARNING: Found \$overlap_count overlapping samples between input and reference" >> ${meta.id}_sample_check.txt
        echo "Overlapping samples:" >> ${meta.id}_sample_check.txt
        head -10 overlapping_samples.txt >> ${meta.id}_sample_check.txt
        if [ \$overlap_count -gt 10 ]; then
            echo "... and \$((overlap_count - 10)) more" >> ${meta.id}_sample_check.txt
        fi
    fi
    echo "  Report generation completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_sample_timing.log
    
    # Final timing report
    END_TIME=\$(date +%s)
    echo "Process completed at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${meta.id}_sample_timing.log
    echo "Total execution time: \$((END_TIME - START_TIME)) seconds" >> ${meta.id}_sample_timing.log
    echo "=== End of Timing Report ===" >> ${meta.id}_sample_timing.log
    
    # Display timing summary
    cat ${meta.id}_sample_timing.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}