process FILTER_MIN_AC {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.minac.vcf.gz"), path("*.minac.vcf.gz.tbi"), emit: vcf
    path "*.filtered.log"                                               , emit: log
    path "versions.yml"                                                 , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_ac = params.min_ac ?: 2
    """
    # Filter by minimum allele count
    bcftools filter \\
        --include "INFO/AC>=${min_ac}" \\
        --output-type z \\
        --output ${prefix}.minac.vcf.gz \\
        ${vcf} \\
        2> ${prefix}.filtered.log
    
    # Index the output
    bcftools index -t ${prefix}.minac.vcf.gz
    
    # Report statistics
    echo "Variants before filtering:"
    bcftools view -H ${vcf} | wc -l
    echo "Variants after filtering (AC>=${min_ac}):"
    bcftools view -H ${prefix}.minac.vcf.gz | wc -l
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.minac.vcf.gz
    touch ${prefix}.minac.vcf.gz.tbi
    touch ${prefix}.filtered.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}