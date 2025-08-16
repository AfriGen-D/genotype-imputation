process CREATE_CHUNK_LIST {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val chunk_size  // Size in base pairs

    output:
    tuple val(meta), path("chunk_list.txt"), emit: chunk_list
    path "chunk_summary.txt", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Creating genomic chunks of ${chunk_size} bp..."
    
    # Get chromosome sizes from the index
    bcftools index -s ${vcf} > chr_stats.txt
    
    # Create genomic chunks based on fixed window size
    echo "chunk_name\tregion" > chunk_list.txt
    
    chunk_id=0
    while IFS=\$'\\t' read -r chr length n_records; do
        if [ "\$n_records" -eq 0 ]; then
            continue
        fi
        
        echo "Processing chromosome \$chr (length: \$length bp, variants: \$n_records)"
        
        # Create chunks for this chromosome
        start=1
        while [ \$start -le \$length ]; do
            end=\$((start + ${chunk_size} - 1))
            if [ \$end -gt \$length ]; then
                end=\$length
            fi
            
            # Create chunk (we'll check for variants during extraction)
            printf "%s_chunk_%04d\\t%s:%d-%d\\n" "${prefix}" \$chunk_id "\$chr" \$start \$end >> chunk_list.txt
            chunk_id=\$((chunk_id + 1))
            
            start=\$((end + 1))
        done
    done < chr_stats.txt
    
    # Create summary
    echo "Chunk creation plan:" > chunk_summary.txt
    total_variants=\$(bcftools index -n ${vcf})
    echo "Total variants: \$total_variants" >> chunk_summary.txt
    echo "Target chunk size: ${chunk_size} bp" >> chunk_summary.txt
    n_chunks=\$(tail -n +2 chunk_list.txt | wc -l)
    echo "Number of chunks: \$n_chunks" >> chunk_summary.txt
    echo "" >> chunk_summary.txt
    echo "Chunk details:" >> chunk_summary.txt
    tail -n +2 chunk_list.txt | head -20 >> chunk_summary.txt
    if [ \$n_chunks -gt 20 ]; then
        echo "... and \$((n_chunks - 20)) more chunks" >> chunk_summary.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}