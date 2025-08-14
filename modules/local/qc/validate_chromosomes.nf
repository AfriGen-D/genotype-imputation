process VALIDATE_CHROMOSOMES {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf)
    val chromosomes
    
    output:
    tuple val(meta), path("${meta.id}_chr_validation.txt"), emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    # Get list of chromosomes in VCF
    bcftools query -f '%CHROM\\n' $vcf | sort -u > vcf_chromosomes.txt
    
    # Generate validation report
    echo "Chromosome Validation for ${meta.id}" > ${meta.id}_chr_validation.txt
    echo "Requested chromosomes: $chromosomes" >> ${meta.id}_chr_validation.txt
    echo "VCF contains chromosomes:" >> ${meta.id}_chr_validation.txt
    cat vcf_chromosomes.txt >> ${meta.id}_chr_validation.txt
    
    # Check if requested chromosomes are present
    if [ "$chromosomes" != "ALL" ]; then
        for chr in \$(echo $chromosomes | tr ',' ' '); do
            if grep -q "^\$chr\$" vcf_chromosomes.txt; then
                echo "Chromosome \$chr: FOUND" >> ${meta.id}_chr_validation.txt
            else
                echo "WARNING: Chromosome \$chr: NOT FOUND" >> ${meta.id}_chr_validation.txt
            fi
        done
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}