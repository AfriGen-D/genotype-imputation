process REDISTRIBUTE_FAILED_CHUNKS {
    tag "$meta.sample"
    label 'process_medium'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), path(vcf_chunks), path(status_files)
    
    output:
    tuple val(meta), path("*.redistributed.vcf.gz"), path("*.redistributed.vcf.gz.tbi"), emit: vcf
    path "*.redistribution.log"                                                        , emit: log
    path "versions.yml"                                                               , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.sample}"
    """
    # Create a redistribution plan
    echo "Analyzing chunks for dataset ${meta.sample}" > ${prefix}.redistribution.log
    
    # Parse status files to identify passed and failed chunks
    touch passed_chunks.txt failed_chunks.txt
    
    for status_file in *.overlap.status; do
        chunk_name=\$(basename \$status_file .overlap.status)
        status=\$(cat \$status_file)
        if [ "\$status" == "PASS" ]; then
            echo "\$chunk_name" >> passed_chunks.txt
        else
            echo "\$chunk_name" >> failed_chunks.txt
        fi
    done
    
    num_passed=\$(wc -l < passed_chunks.txt)
    num_failed=\$(wc -l < failed_chunks.txt)
    
    echo "Passed chunks: \$num_passed" >> ${prefix}.redistribution.log
    echo "Failed chunks: \$num_failed" >> ${prefix}.redistribution.log
    
    if [ \$num_failed -eq 0 ]; then
        echo "No failed chunks to redistribute" >> ${prefix}.redistribution.log
        # Just combine all passed chunks
        bcftools concat \\
            --allow-overlaps \\
            --output-type z \\
            --output ${prefix}.redistributed.vcf.gz \\
            --threads ${task.cpus} \\
            \$(cat passed_chunks.txt | sed 's/\$/.vcf.gz/g')
    elif [ \$num_passed -eq 0 ]; then
        echo "ERROR: All chunks failed overlap check!" >> ${prefix}.redistribution.log
        echo "Cannot proceed with redistribution" >> ${prefix}.redistribution.log
        exit 1
    else
        echo "Redistributing \$num_failed failed chunks to \$num_passed passing chunks" >> ${prefix}.redistribution.log
        
        # For each failed chunk, find the nearest passing chunk and merge
        while read failed_chunk; do
            # Extract chromosome and position from chunk name
            # Format: dataset_chr_start_end
            chr=\$(echo \$failed_chunk | cut -d'_' -f2)
            start=\$(echo \$failed_chunk | cut -d'_' -f3)
            end=\$(echo \$failed_chunk | cut -d'_' -f4)
            
            echo "Processing failed chunk: \$failed_chunk (chr\$chr:\$start-\$end)" >> ${prefix}.redistribution.log
            
            # Find the nearest passing chunk on the same chromosome
            nearest_chunk=""
            min_distance=999999999999
            
            while read passed_chunk; do
                p_chr=\$(echo \$passed_chunk | cut -d'_' -f2)
                p_start=\$(echo \$passed_chunk | cut -d'_' -f3)
                p_end=\$(echo \$passed_chunk | cut -d'_' -f4)
                
                # Only consider chunks on the same chromosome
                if [ "\$chr" == "\$p_chr" ]; then
                    # Calculate distance (use start positions)
                    if [ \$start -ge \$p_start ]; then
                        distance=\$((start - p_start))
                    else
                        distance=\$((p_start - start))
                    fi
                    
                    if [ \$distance -lt \$min_distance ]; then
                        min_distance=\$distance
                        nearest_chunk=\$passed_chunk
                    fi
                fi
            done < passed_chunks.txt
            
            if [ -n "\$nearest_chunk" ]; then
                echo "  -> Merging with nearest chunk: \$nearest_chunk (distance: \$min_distance bp)" >> ${prefix}.redistribution.log
                
                # Merge the failed chunk variants with the nearest passing chunk
                bcftools concat \\
                    --allow-overlaps \\
                    --remove-duplicates \\
                    --output-type z \\
                    --output \${nearest_chunk}.merged.vcf.gz \\
                    \${nearest_chunk}.vcf.gz \${failed_chunk}.vcf.gz
                
                # Replace the original passing chunk with the merged version
                mv \${nearest_chunk}.merged.vcf.gz \${nearest_chunk}.vcf.gz
                bcftools index -t \${nearest_chunk}.vcf.gz
            else
                echo "  -> WARNING: No passing chunk found on chromosome \$chr for failed chunk \$failed_chunk" >> ${prefix}.redistribution.log
            fi
        done < failed_chunks.txt
        
        # Combine all redistributed chunks
        echo "Creating final redistributed VCF..." >> ${prefix}.redistribution.log
        bcftools concat \\
            --allow-overlaps \\
            --remove-duplicates \\
            --output-type z \\
            --output ${prefix}.redistributed.vcf.gz \\
            --threads ${task.cpus} \\
            \$(cat passed_chunks.txt | sed 's/\$/.vcf.gz/g')
    fi
    
    # Index the output
    bcftools index -t ${prefix}.redistributed.vcf.gz
    
    # Report final statistics
    echo "" >> ${prefix}.redistribution.log
    echo "Final statistics:" >> ${prefix}.redistribution.log
    echo "Total variants in redistributed VCF: \$(bcftools view -H ${prefix}.redistributed.vcf.gz | wc -l)" >> ${prefix}.redistribution.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.sample}"
    """
    touch ${prefix}.redistributed.vcf.gz
    touch ${prefix}.redistributed.vcf.gz.tbi
    touch ${prefix}.redistribution.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}