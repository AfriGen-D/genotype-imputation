process REMOVE_DUPLICATES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.no_dups.vcf.gz"), path("${prefix}.no_dups.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.duplicates.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract variant IDs and check for duplicates
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${vcf} | \\
        sort | uniq -d > ${prefix}.duplicates.txt
    
    # Count duplicates
    n_dups=\$(wc -l < ${prefix}.duplicates.txt)
    echo "Found \$n_dups duplicate variants" | tee -a ${prefix}.duplicates.txt
    
    # Remove duplicates using bcftools norm
    bcftools norm \\
        --rm-dup all \\
        --output-type z \\
        --output ${prefix}.no_dups.vcf.gz \\
        --threads ${task.cpus} \\
        ${args} \\
        ${vcf}
    
    # Index output
    bcftools index -t ${prefix}.no_dups.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}