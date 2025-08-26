process EXTRACT_VCF_CHUNK {
    tag "$chunk_name"
    label 'process_low'

    input:
    tuple val(meta), val(chunk_name), val(region), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("${chunk_name}.vcf.gz"), path("${chunk_name}.vcf.gz.tbi"), emit: chunk
    path "${chunk_name}.stats.txt", emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Extracting chunk ${chunk_name} for region ${region}"
    
    # Extract chunk using the region
    bcftools view \\
        -r ${region} \\
        ${vcf} \\
        -Oz \\
        -o ${chunk_name}.vcf.gz \\
        --threads ${task.cpus}
    
    # Check if chunk has variants
    n_vars=\$(bcftools view -H ${chunk_name}.vcf.gz | wc -l)
    
    if [ \$n_vars -eq 0 ]; then
        echo "WARNING: Chunk ${chunk_name} has no variants in region ${region}"
        # Create a minimal valid VCF with just headers
        bcftools view -h ${vcf} | bgzip -c > ${chunk_name}.vcf.gz
    fi
    
    # Index the chunk
    bcftools index -t ${chunk_name}.vcf.gz --threads ${task.cpus}
    
    # Get statistics
    n_vars=\$(bcftools index -n ${chunk_name}.vcf.gz)
    echo "Chunk: ${chunk_name}" > ${chunk_name}.stats.txt
    echo "Region: ${region}" >> ${chunk_name}.stats.txt
    echo "Variants: \$n_vars" >> ${chunk_name}.stats.txt
    
    echo "Created ${chunk_name}.vcf.gz with \$n_vars variants"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}