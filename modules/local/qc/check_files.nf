process CHECK_FILES {
    tag "$meta.id"
    label 'process_single'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path(vcf), emit: vcf
    path "versions.yml"        , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "=== Validating required files for sample ${meta.id} ==="
    
    # Check if VCF file exists and is valid
    if [ ! -f "${vcf}" ]; then
        echo "ERROR: VCF file ${vcf} not found"
        exit 1
    fi
    
    # Check VCF integrity
    bcftools view -h ${vcf} > /dev/null 2>&1
    if [ \$? -ne 0 ]; then
        echo "ERROR: ${vcf} is not a valid VCF file"
        exit 1
    fi
    
    # Check if VCF has any variants (handle SIGPIPE from head)
    set +e  # Temporarily allow non-zero exit codes
    variant_count=\$(bcftools view -H ${vcf} 2>/dev/null | head -1 | wc -l)
    pipe_status=\${PIPESTATUS[0]}
    set -e  # Re-enable exit on error
    
    # Check if bcftools failed for reasons other than SIGPIPE (141)
    if [ "\$pipe_status" -ne 0 ] && [ "\$pipe_status" -ne 141 ]; then
        echo "ERROR: Failed to read variants from ${vcf}"
        exit 1
    fi
    
    if [ "\$variant_count" -eq 0 ]; then
        echo "ERROR: VCF file ${vcf} is empty (contains no variants)"
        echo "       The pipeline requires VCF files with actual variant data to process."
        echo "       Please provide a VCF file with variant records."
        exit 1
    fi
    
    echo "✓ VCF file ${vcf} validated successfully"
    
    # Check Eagle genetic map file (if provided via params)
    if [ ! -z "${params.eagle_genetic_map}" ] && [ "${params.eagle_genetic_map}" != "null" ]; then
        if [ ! -f "${params.eagle_genetic_map}" ]; then
            echo "ERROR: Eagle genetic map file not found: ${params.eagle_genetic_map}"
            echo "Please ensure the file exists or update the eagle_genetic_map parameter"
            exit 1
        fi
        echo "✓ Eagle genetic map file found: ${params.eagle_genetic_map}"
    else
        echo "WARNING: Eagle genetic map not provided - will be required for phasing"
    fi
    
    # Check reference genome (if provided)
    if [ ! -z "${params.reference_genome}" ] && [ "${params.reference_genome}" != "null" ]; then
        if [ ! -f "${params.reference_genome}" ]; then
            echo "ERROR: Reference genome file not found: ${params.reference_genome}"
            echo "Please ensure the file exists or update the reference_genome parameter"
            exit 1
        fi
        echo "✓ Reference genome file found: ${params.reference_genome}"
    else
        echo "WARNING: Reference genome not provided - will be required for phasing"
    fi
    
    # Check reference panels (if provided)
    if [ ! -z "${params.ref_panels}" ] && [ "${params.ref_panels}" != "[]" ]; then
        echo "✓ Reference panels configured: ${params.ref_panels}"
        
        # For the nf-core workflow, we need actual reference panel VCF files
        # These should be passed properly through the workflow
        # For now, check if the reference panel paths exist using the first panel as example
        # TODO: Properly pass reference panel VCFs to phasing and imputation steps
        echo "INFO: Reference panels need to be properly integrated into phasing workflow"
    else
        echo "WARNING: No reference panels configured - will be required for imputation"
    fi
    
    echo "=== File validation complete for ${meta.id} ==="
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Stub: Checking ${vcf}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}