process FILTER_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    val min_ac
    val min_maf
    val max_missing

    output:
    tuple val(meta), path("${prefix}.filtered.vcf.gz"), path("${prefix}.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.filter_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Count variants before filtering
    n_before=\$(bcftools view -H ${vcf} | wc -l)
    echo "Variants before filtering: \$n_before" > ${prefix}.filter_stats.txt
    
    # Apply filters
    bcftools view \\
        --min-ac ${min_ac} \\
        --min-af ${min_maf} \\
        --max-missing ${max_missing} \\
        --output-type z \\
        --output ${prefix}.filtered.vcf.gz \\
        --threads ${task.cpus} \\
        ${args} \\
        ${vcf}
    
    # Index output
    bcftools index -t ${prefix}.filtered.vcf.gz
    
    # Count variants after filtering
    n_after=\$(bcftools view -H ${prefix}.filtered.vcf.gz | wc -l)
    echo "Variants after filtering: \$n_after" >> ${prefix}.filter_stats.txt
    echo "Variants removed: \$((n_before - n_after))" >> ${prefix}.filter_stats.txt
    echo "Removal rate: \$(echo "scale=2; 100 * (n_before - n_after) / n_before" | bc)%" >> ${prefix}.filter_stats.txt
    
    # Get MAF distribution
    bcftools +fill-tags ${prefix}.filtered.vcf.gz -- -t AF | \\
        bcftools query -f '%AF\\n' | \\
        awk '{
            if(\$1 < 0.01) maf_bins["0-0.01"]++
            else if(\$1 < 0.05) maf_bins["0.01-0.05"]++
            else if(\$1 < 0.1) maf_bins["0.05-0.1"]++
            else if(\$1 < 0.5) maf_bins["0.1-0.5"]++
            else maf_bins["0.5-1"]++
        } END {
            print "\\nMAF distribution:"
            for(bin in maf_bins) print bin": "maf_bins[bin]
        }' >> ${prefix}.filter_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}