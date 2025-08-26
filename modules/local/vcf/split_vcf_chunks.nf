process SPLIT_VCF_CHUNKS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val chunk_size

    output:
    tuple val(meta), path("${meta.id}_chunk_*.vcf.gz"), path("${meta.id}_chunk_*.vcf.gz.tbi"), emit: chunks
    path "chunk_regions.txt", emit: regions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Creating chunks of ${chunk_size} variants each..."
    
    # Get total variant count
    total_variants=\$(bcftools index -n ${vcf})
    echo "Total variants: \$total_variants"
    
    # Calculate number of chunks
    num_chunks=\$(( (total_variants + ${chunk_size} - 1) / ${chunk_size} ))
    echo "Creating \$num_chunks chunks of up to ${chunk_size} variants each"
    
    # Get chromosome info
    bcftools index -s ${vcf} > chr_stats.txt
    
    # Simple approach: split VCF into equal-sized chunks by variant count
    # Using awk to create chunk regions
    
    chunk_id=0
    current_count=0
    target_per_chunk=${chunk_size}
    
    # Create chunks based on variant count
    bcftools query -f '%CHROM\\t%POS\\n' ${vcf} | \\
    awk -v chunk_size=${chunk_size} -v prefix="${prefix}" '
    BEGIN {
        chunk_id = 0
        chunk_start = ""
        chunk_chr = ""
        count = 0
    }
    {
        if (count == 0) {
            # Start new chunk
            chunk_start = \$2
            chunk_chr = \$1
        }
        
        count++
        
        if (count >= chunk_size || \$1 != chunk_chr) {
            # End current chunk
            if (chunk_chr != "") {
                printf "%s_chunk_%04d\\t%s:%s-%s\\n", prefix, chunk_id, chunk_chr, chunk_start, prev_pos
                chunk_id++
            }
            
            # Start new chunk if chromosome changed
            if (\$1 != chunk_chr) {
                chunk_start = \$2
                chunk_chr = \$1
                count = 1
            } else {
                count = 0
            }
        }
        
        prev_pos = \$2
        prev_chr = \$1
    }
    END {
        # Output last chunk
        if (count > 0) {
            printf "%s_chunk_%04d\\t%s:%s-%s\\n", prefix, chunk_id, prev_chr, chunk_start, prev_pos
        }
    }
    ' > chunk_list.txt
    
    # Create VCF chunks from the list
    touch chunk_regions.txt
    while IFS=\$'\\t' read -r chunk_name region; do
        echo "Creating \$chunk_name for region \$region"
        echo "\$chunk_name\\t\$region" >> chunk_regions.txt
        
        # Extract chunk
        bcftools view -r \$region ${vcf} -Oz -o \${chunk_name}.vcf.gz
        
        # Index chunk
        bcftools index -t \${chunk_name}.vcf.gz
        
        # Count variants
        n_vars=\$(bcftools index -n \${chunk_name}.vcf.gz)
        echo "  Created \${chunk_name}.vcf.gz with \$n_vars variants"
    done < chunk_list.txt
    
    # Verify we created chunks
    n_chunks=\$(ls ${prefix}_chunk_*.vcf.gz 2>/dev/null | wc -l)
    if [ "\$n_chunks" -eq 0 ]; then
        echo "Warning: No chunks created, using original file as single chunk"
        cp ${vcf} ${prefix}_chunk_0000.vcf.gz
        cp ${vcf_index} ${prefix}_chunk_0000.vcf.gz.tbi
        echo "${prefix}_chunk_0000\\tfull" > chunk_regions.txt
    fi
    
    echo "Successfully created \$n_chunks chunks"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}