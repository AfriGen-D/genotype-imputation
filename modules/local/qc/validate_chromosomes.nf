process VALIDATE_CHROMOSOMES {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf)
    val chromosomes
    
    output:
    tuple val(meta), path("${meta.id}_chr_validation.txt"), emit: report
    tuple val(meta), path("${meta.id}_chr_timing.log"), emit: timing
    path "versions.yml", emit: versions
    
    script:
    """
    # Initialize timing log
    echo "=== VALIDATE_CHROMOSOMES Timing Report ===" > ${meta.id}_chr_timing.log
    echo "Process started at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${meta.id}_chr_timing.log
    START_TIME=\$(date +%s)
    
    # Get list of chromosomes in VCF
    echo "Step 1: Extracting chromosomes from VCF started at \$(date '+%H:%M:%S')" >> ${meta.id}_chr_timing.log
    STEP_START=\$(date +%s)
    bcftools query -f '%CHROM\\n' $vcf | sort -u > vcf_chromosomes.txt
    echo "  Chromosome extraction completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_chr_timing.log
    
    # Generate validation report
    echo "Step 2: Generating validation report started at \$(date '+%H:%M:%S')" >> ${meta.id}_chr_timing.log
    STEP_START=\$(date +%s)
    echo "Chromosome Validation for ${meta.id}" > ${meta.id}_chr_validation.txt
    echo "Requested chromosomes: $chromosomes" >> ${meta.id}_chr_validation.txt
    echo "VCF contains chromosomes:" >> ${meta.id}_chr_validation.txt
    cat vcf_chromosomes.txt >> ${meta.id}_chr_validation.txt
    echo "  Report generation completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_chr_timing.log
    
    # Check if requested chromosomes are present
    echo "Step 3: Chromosome validation started at \$(date '+%H:%M:%S')" >> ${meta.id}_chr_timing.log
    STEP_START=\$(date +%s)
    if [ "$chromosomes" != "ALL" ]; then
        for chr in \$(echo $chromosomes | tr ',' ' '); do
            if grep -q "^\$chr\$" vcf_chromosomes.txt; then
                echo "Chromosome \$chr: FOUND" >> ${meta.id}_chr_validation.txt
            else
                echo "WARNING: Chromosome \$chr: NOT FOUND" >> ${meta.id}_chr_validation.txt
            fi
        done
    fi
    echo "  Chromosome validation completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${meta.id}_chr_timing.log
    
    # Final timing report
    END_TIME=\$(date +%s)
    echo "Process completed at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${meta.id}_chr_timing.log
    echo "Total execution time: \$((END_TIME - START_TIME)) seconds" >> ${meta.id}_chr_timing.log
    echo "=== End of Timing Report ===" >> ${meta.id}_chr_timing.log
    
    # Display timing summary
    cat ${meta.id}_chr_timing.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}