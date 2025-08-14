process FILTER_IMPUTED {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    val info_cutoff
    val maf_cutoff

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: filtered_vcf
    tuple val(meta), path("*.filter_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def info_filter = info_cutoff > 0 ? "INFO>=${info_cutoff}" : ''
    def maf_filter = maf_cutoff > 0 ? "MAF>=${maf_cutoff}" : ''
    def filters = [info_filter, maf_filter].findAll { it }.join(' && ')
    def filter_expression = filters ? "--include '${filters}'" : ''
    """
    # Count variants before filtering
    TOTAL_BEFORE=\$(bcftools view -H ${vcf} | wc -l)

    # Apply filters
    bcftools view \\
        ${filter_expression} \\
        -O z \\
        -o ${prefix}.filtered.vcf.gz \\
        ${vcf} \\
        ${args}

    # Index filtered VCF
    bcftools index -t ${prefix}.filtered.vcf.gz

    # Count variants after filtering
    TOTAL_AFTER=\$(bcftools view -H ${prefix}.filtered.vcf.gz | wc -l)

    # Generate filter statistics
    echo "# Filter Statistics for ${prefix}" > ${prefix}.filter_stats.txt
    echo "# INFO cutoff: ${info_cutoff}" >> ${prefix}.filter_stats.txt
    echo "# MAF cutoff: ${maf_cutoff}" >> ${prefix}.filter_stats.txt
    echo "" >> ${prefix}.filter_stats.txt
    echo "Variants before filtering: \$TOTAL_BEFORE" >> ${prefix}.filter_stats.txt
    echo "Variants after filtering: \$TOTAL_AFTER" >> ${prefix}.filter_stats.txt
    
    if [ \$TOTAL_BEFORE -gt 0 ]; then
        RETAINED_PCT=\$(echo "scale=2; (\$TOTAL_AFTER / \$TOTAL_BEFORE) * 100" | bc -l)
        REMOVED=\$((\$TOTAL_BEFORE - \$TOTAL_AFTER))
        REMOVED_PCT=\$(echo "scale=2; (\$REMOVED / \$TOTAL_BEFORE) * 100" | bc -l)
        echo "Variants retained: \$TOTAL_AFTER (\${RETAINED_PCT}%)" >> ${prefix}.filter_stats.txt
        echo "Variants removed: \$REMOVED (\${REMOVED_PCT}%)" >> ${prefix}.filter_stats.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi
    touch ${prefix}.filter_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}