process SITES_ONLY {
    tag "$meta.id"
    label 'process_single'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.sites.vcf.gz"), path("*.sites.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.sites.txt")                                , emit: sites
    path "versions.yml"                                                 , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract sites only (remove genotype information)
    bcftools view \\
        --drop-genotypes \\
        --output-type z \\
        --output ${prefix}.sites.vcf.gz \\
        ${vcf}
    
    # Index the output
    bcftools index -t ${prefix}.sites.vcf.gz
    
    # Create sites list
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' \\
        ${prefix}.sites.vcf.gz > ${prefix}.sites.txt
    
    # Report number of sites
    echo "Total sites: \$(wc -l < ${prefix}.sites.txt)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sites.vcf.gz
    touch ${prefix}.sites.vcf.gz.tbi
    touch ${prefix}.sites.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}