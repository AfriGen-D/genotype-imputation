process GENERATE_CHUNK_MAP {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), path(chunk_vcf)
    
    output:
    tuple val(meta), path("*.chunk.map"), emit: map
    path "versions.yml"                  , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Generate map file with chromosome and position for chunk
    bcftools query -f '%CHROM\\t%POS\\n' ${chunk_vcf} > ${prefix}.chunk.map
    
    # Report statistics
    echo "Generated map file for chunk ${prefix}"
    echo "Total variants: \$(wc -l < ${prefix}.chunk.map)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chunk.map
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}