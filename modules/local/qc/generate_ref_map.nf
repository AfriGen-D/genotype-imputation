process GENERATE_REF_MAP {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), path(reference_panel), path(reference_panel_index)
    
    output:
    tuple val(meta), path("*.ref.map"), emit: map
    path "versions.yml"                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chrm = meta.contig ?: 'chr21'
    def chunk_start = meta.start ?: ''
    def chunk_end = meta.end ?: ''
    """
    # Generate map file with chromosome and position for reference panel in chunk region
    if [ -n "${chunk_start}" ] && [ -n "${chunk_end}" ]; then
        bcftools query -r ${chrm}:${chunk_start}-${chunk_end} -f '%CHROM\\t%POS\\n' ${reference_panel} > ${prefix}.ref.map
    else
        bcftools query -f '%CHROM\\t%POS\\n' ${reference_panel} > ${prefix}.ref.map
    fi
    
    # Report statistics
    echo "Generated map file for reference panel ${prefix}"
    echo "Region: ${chrm}:${chunk_start}-${chunk_end}"
    echo "Total variants: \$(wc -l < ${prefix}.ref.map)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ref.map
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}