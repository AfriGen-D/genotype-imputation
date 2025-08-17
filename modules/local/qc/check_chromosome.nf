process CHECK_CHROMOSOME {
    tag "$meta.id"
    label 'process_single'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path(vcf)         , emit: vcf
    tuple val(meta), path("*.chroms")  , emit: chromosomes
    path "versions.yml"                 , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create index if it doesn't exist
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        bcftools index -t ${vcf}
    fi
    
    # Extract chromosome list from VCF
    bcftools index -s ${vcf} | cut -f1 > ${prefix}.chroms
    
    # Check if chromosomes match expected format
    echo "Chromosomes found in ${vcf}:"
    cat ${prefix}.chroms
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chroms
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}