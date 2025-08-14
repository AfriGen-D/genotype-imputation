process CHECK_REFERENCE_COMPATIBILITY {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf)
    tuple val(ref_name), path(ref_vcf), path(ref_m3vcf)
    
    output:
    tuple val(meta), path("${meta.id}_ref_check.txt"), emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    # Check chromosome naming convention
    bcftools query -f '%CHROM\\n' $vcf | head -1 > chr_format.txt
    bcftools query -f '%CHROM\\n' $ref_vcf | head -1 > ref_chr_format.txt
    
    # Generate compatibility report
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
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}