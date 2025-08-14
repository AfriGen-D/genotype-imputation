process MERGE_IMPUTED_CHUNKS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: merged_vcf
    tuple val(meta), path("*.merge_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create file list for merging
    ls *.vcf.gz | sort -V > vcf_list.txt

    # Log merge information
    echo "# Merge log for ${prefix}" > ${prefix}.merge_log.txt
    echo "# Date: \$(date)" >> ${prefix}.merge_log.txt
    echo "# Files to merge:" >> ${prefix}.merge_log.txt
    cat vcf_list.txt >> ${prefix}.merge_log.txt
    echo "" >> ${prefix}.merge_log.txt

    # Count variants in each file
    echo "# Variant counts per file:" >> ${prefix}.merge_log.txt
    while read vcf; do
        count=\$(bcftools view -H \$vcf | wc -l)
        echo "\$vcf: \$count variants" >> ${prefix}.merge_log.txt
    done < vcf_list.txt
    echo "" >> ${prefix}.merge_log.txt

    # Merge VCF files
    bcftools concat \\
        --file-list vcf_list.txt \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        --threads ${task.cpus} \\
        ${args}

    # Index merged VCF
    bcftools index -t ${prefix}.merged.vcf.gz

    # Final statistics
    TOTAL_MERGED=\$(bcftools view -H ${prefix}.merged.vcf.gz | wc -l)
    echo "# Final statistics:" >> ${prefix}.merge_log.txt
    echo "Total variants in merged file: \$TOTAL_MERGED" >> ${prefix}.merge_log.txt
    
    # Check for chromosome coverage
    echo "# Chromosome coverage:" >> ${prefix}.merge_log.txt
    bcftools index -s ${prefix}.merged.vcf.gz >> ${prefix}.merge_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi
    touch ${prefix}.merge_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}