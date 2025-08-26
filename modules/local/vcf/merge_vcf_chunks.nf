process MERGE_VCF_CHUNKS {
    tag "$meta.id"
    label 'process_medium'
    
    input:
    tuple val(meta), path(chunk_vcfs)
    
    output:
    tuple val(meta), path("${meta.id}.merged.vcf.gz"), path("${meta.id}.merged.vcf.gz.tbi"), emit: vcf
    path "merge_stats.txt", emit: stats
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Merging ${chunk_vcfs.size()} VCF chunks..."
    
    # Sort chunks by name to ensure correct order
    ls ${meta.id}_chunk_*.vcf.gz | sort > chunks_list.txt
    
    # Count total variants before merge
    TOTAL_BEFORE=0
    for vcf in \$(cat chunks_list.txt); do
        COUNT=\$(bcftools query -f 'x\\n' \$vcf | wc -l)
        echo "\$vcf: \$COUNT variants" >> merge_stats.txt
        TOTAL_BEFORE=\$((TOTAL_BEFORE + COUNT))
    done
    echo "Total variants before merge: \$TOTAL_BEFORE" >> merge_stats.txt
    
    # Merge all chunks
    if [ \$(cat chunks_list.txt | wc -l) -eq 1 ]; then
        # Only one chunk, just rename
        cp \$(cat chunks_list.txt) ${prefix}.merged.vcf.gz
    else
        # Multiple chunks, concatenate
        bcftools concat \\
            --file-list chunks_list.txt \\
            --allow-overlaps \\
            --remove-duplicates \\
            -Oz -o ${prefix}.merged.vcf.gz \\
            ${args}
    fi
    
    # Index merged file
    bcftools index -t ${prefix}.merged.vcf.gz
    
    # Count variants after merge
    TOTAL_AFTER=\$(bcftools query -f 'x\\n' ${prefix}.merged.vcf.gz | wc -l)
    echo "Total variants after merge: \$TOTAL_AFTER" >> merge_stats.txt
    
    # Verify merge integrity
    if [ \$TOTAL_AFTER -ne \$TOTAL_BEFORE ]; then
        echo "WARNING: Variant count mismatch (before: \$TOTAL_BEFORE, after: \$TOTAL_AFTER)" >> merge_stats.txt
        echo "This may be due to duplicate removal" >> merge_stats.txt
    fi
    
    echo "Merge completed successfully"
    cat merge_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}