process EXTRACT_INFO_SCORES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)
    val info_cutoff

    output:
    tuple val(meta), path("*.info.gz"), emit: info
    tuple val(meta), path("*.summary.txt"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract info scores from VCF
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/INFO\\t%INFO/MAF\\n' \\
        ${vcf} | gzip > ${prefix}.info.gz

    # Generate summary statistics
    echo "# INFO Score Summary for ${prefix}" > ${prefix}.summary.txt
    echo "# Cutoff: ${info_cutoff}" >> ${prefix}.summary.txt
    echo "" >> ${prefix}.summary.txt

    zcat ${prefix}.info.gz | awk -v cutoff=${info_cutoff} '
    BEGIN {
        total = 0
        well_imputed = 0
        sum_info = 0
        min_info = 1
        max_info = 0
    }
    NR > 1 && \$6 != "." {
        total++
        info = \$6
        sum_info += info
        
        if (info >= cutoff) well_imputed++
        if (info < min_info) min_info = info
        if (info > max_info) max_info = info
    }
    END {
        avg_info = (total > 0) ? sum_info / total : 0
        pct_well_imputed = (total > 0) ? (well_imputed / total) * 100 : 0
        
        print "Total variants: " total
        print "Well-imputed variants (INFO >= " cutoff "): " well_imputed " (" pct_well_imputed "%)"
        print "Average INFO score: " avg_info
        print "Min INFO score: " min_info
        print "Max INFO score: " max_info
    }' >> ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.info.gz
    touch ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}