process QC_SITE_MISSINGNESS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.nomiss.vcf.gz"), path("*.nomiss.vcf.gz.tbi"), emit: vcf
    path "*.missingness.txt"                                              , emit: missingness
    path "versions.yml"                                                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_missing = params.site_missingness ?: 0.05
    """
    # Calculate missingness per site
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AN\\t%INFO/AC\\n' ${vcf} > ${prefix}.info.txt
    
    # Filter sites with high missingness
    bcftools filter \\
        --include "F_MISSING<${max_missing}" \\
        --output-type z \\
        --output ${prefix}.nomiss.vcf.gz \\
        ${vcf}
    
    # Index the output
    bcftools index -t ${prefix}.nomiss.vcf.gz
    
    # Report missingness statistics
    echo "Site missingness statistics:" > ${prefix}.missingness.txt
    echo "Threshold: ${max_missing}" >> ${prefix}.missingness.txt
    echo "Sites before filtering:" >> ${prefix}.missingness.txt
    bcftools view -H ${vcf} | wc -l >> ${prefix}.missingness.txt
    echo "Sites after filtering:" >> ${prefix}.missingness.txt
    bcftools view -H ${prefix}.nomiss.vcf.gz | wc -l >> ${prefix}.missingness.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nomiss.vcf.gz
    touch ${prefix}.nomiss.vcf.gz.tbi
    touch ${prefix}.missingness.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}