process CHECK_SAMPLE_OVERLAP {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf)
    tuple val(ref_name), path(ref_vcf), path(ref_m3vcf)
    
    output:
    tuple val(meta), path("${meta.id}_sample_check.txt"), emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    # Get sample lists
    bcftools query -l $vcf > vcf_samples.txt
    bcftools query -l $ref_vcf > ref_samples.txt 2>/dev/null || touch ref_samples.txt
    
    # Count samples
    vcf_count=\$(wc -l < vcf_samples.txt)
    ref_count=\$(wc -l < ref_samples.txt)
    
    # Check for overlapping samples
    comm -12 <(sort vcf_samples.txt) <(sort ref_samples.txt) > overlapping_samples.txt || true
    overlap_count=\$(wc -l < overlapping_samples.txt)
    
    # Generate report
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
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}