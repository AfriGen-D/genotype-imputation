process QC_DUPL {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.nodup.vcf.gz"), path("*.nodup.vcf.gz.tbi"), emit: vcf
    path "*.duplicates.txt"                                            , emit: duplicates
    path "versions.yml"                                                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Remove duplicate variants
    bcftools norm \\
        --rm-dup all \\
        --output-type z \\
        --output ${prefix}.nodup.vcf.gz \\
        ${vcf} \\
        2> ${prefix}.duplicates.txt
    
    # Index the output
    bcftools index -t ${prefix}.nodup.vcf.gz
    
    # Report duplicates found
    echo "Duplicate variants removed:"
    grep -c "duplicate" ${prefix}.duplicates.txt || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nodup.vcf.gz
    touch ${prefix}.nodup.vcf.gz.tbi
    touch ${prefix}.duplicates.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}