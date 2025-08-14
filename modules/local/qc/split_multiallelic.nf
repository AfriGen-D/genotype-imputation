process SPLIT_MULTIALLELIC {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.split.vcf.gz"), path("${prefix}.split.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.split_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Count multiallelic sites before splitting
    n_multi_before=\$(bcftools view -m 3 ${vcf} | grep -v "^#" | wc -l)
    echo "Multiallelic sites before splitting: \$n_multi_before" > ${prefix}.split_stats.txt
    
    # Split multiallelic sites to biallelic
    bcftools norm \\
        --multiallelics -both \\
        --output-type z \\
        --output ${prefix}.split.vcf.gz \\
        --threads ${task.cpus} \\
        ${args} \\
        ${vcf}
    
    # Index output
    bcftools index -t ${prefix}.split.vcf.gz
    
    # Count variants after splitting
    n_vars_after=\$(bcftools view -H ${prefix}.split.vcf.gz | wc -l)
    echo "Total variants after splitting: \$n_vars_after" >> ${prefix}.split_stats.txt
    echo "New variants created: \$((n_vars_after - n_multi_before))" >> ${prefix}.split_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}