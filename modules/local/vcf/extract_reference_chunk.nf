process EXTRACT_REFERENCE_CHUNK {
    tag "$chunk_name:$ref_name"
    label 'process_low'

    input:
    tuple val(chunk_name), val(region), val(ref_name), path(ref_vcf), path(ref_m3vcf)

    output:
    tuple val(chunk_name), val(ref_name), path("${chunk_name}_${ref_name}_ref.vcf.gz"), path("${chunk_name}_${ref_name}_ref.vcf.gz.tbi"), emit: ref_chunk
    path "${chunk_name}_${ref_name}_ref.stats.txt", emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Extracting reference panel chunk for ${chunk_name}, region ${region}, panel ${ref_name}"
    
    # Check if reference VCF has index
    if [ ! -f "${ref_vcf}.tbi" ] && [ ! -f "${ref_vcf}.csi" ]; then
        echo "Creating index for reference panel..."
        bcftools index -t ${ref_vcf} --threads ${task.cpus}
    fi
    
    # Extract reference chunk for the region
    bcftools view \\
        -r ${region} \\
        ${ref_vcf} \\
        -Oz \\
        -o ${chunk_name}_${ref_name}_ref.vcf.gz \\
        --threads ${task.cpus}
    
    # Index the reference chunk
    bcftools index -t ${chunk_name}_${ref_name}_ref.vcf.gz --threads ${task.cpus}
    
    # Get statistics
    n_vars=\$(bcftools index -n ${chunk_name}_${ref_name}_ref.vcf.gz)
    n_samples=\$(bcftools query -l ${chunk_name}_${ref_name}_ref.vcf.gz | wc -l)
    
    echo "Reference chunk: ${chunk_name}" > ${chunk_name}_${ref_name}_ref.stats.txt
    echo "Reference panel: ${ref_name}" >> ${chunk_name}_${ref_name}_ref.stats.txt
    echo "Region: ${region}" >> ${chunk_name}_${ref_name}_ref.stats.txt
    echo "Variants: \$n_vars" >> ${chunk_name}_${ref_name}_ref.stats.txt
    echo "Samples: \$n_samples" >> ${chunk_name}_${ref_name}_ref.stats.txt
    
    echo "Created reference chunk with \$n_vars variants and \$n_samples samples"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}