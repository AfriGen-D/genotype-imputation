process SPLIT_MULTI_ALLELIC {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.split.vcf.gz"), path("*.split.vcf.gz.tbi"), emit: vcf
    path "*.split.log"                                                  , emit: log
    path "versions.yml"                                                 , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Split multi-allelic variants
    bcftools norm \\
        --multiallelics -any \\
        --output-type z \\
        --output ${prefix}.split.vcf.gz \\
        ${vcf} \\
        2> ${prefix}.split.log
    
    # Index the output
    bcftools index -t ${prefix}.split.vcf.gz
    
    # Report statistics
    echo "Multi-allelic variants split:"
    grep -c "split" ${prefix}.split.log || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.split.vcf.gz
    touch ${prefix}.split.vcf.gz.tbi
    touch ${prefix}.split.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}