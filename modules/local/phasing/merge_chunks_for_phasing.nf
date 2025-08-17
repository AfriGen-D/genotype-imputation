process MERGE_CHUNKS_FOR_PHASING {
    tag "$meta.id"
    label 'process_single'
    
    container 'quay.io/biocontainers/bcftools:1.11--h7c999a4_0'
    
    input:
    tuple val(meta), path(vcf1), path(vcf1_index), path(vcf2), path(vcf2_index)
    
    output:
    tuple val(meta), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: merged
    path "versions.yml"                                                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merged_chunk = "${meta.contig}_${meta.start}_${meta.end_extended}"
    """
    # Merge adjacent chunks to create larger region for phasing
    bcftools concat \\
        ${vcf1} ${vcf2} \\
        --allow-overlaps \\
        --remove-duplicates \\
        -Oz -o ${prefix}_${merged_chunk}.merged.vcf.gz
    
    # Index the merged VCF
    bcftools index -t ${prefix}_${merged_chunk}.merged.vcf.gz
    
    # Log merge info
    echo "Merged chunks for low-recombination region:" > ${prefix}_${merged_chunk}.merge_info.txt
    echo "Chunk 1: ${meta.chunk1_id}" >> ${prefix}_${merged_chunk}.merge_info.txt
    echo "Chunk 2: ${meta.chunk2_id}" >> ${prefix}_${merged_chunk}.merge_info.txt
    echo "Merged region: ${meta.contig}:${meta.start}-${meta.end_extended}" >> ${prefix}_${merged_chunk}.merge_info.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merged_chunk = "${meta.contig}_${meta.start}_${meta.end_extended}"
    """
    touch ${prefix}_${merged_chunk}.merged.vcf.gz
    touch ${prefix}_${merged_chunk}.merged.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}