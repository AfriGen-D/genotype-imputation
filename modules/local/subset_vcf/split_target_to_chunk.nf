process SPLIT_TARGET_TO_CHUNK {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(meta_chunk), path(chunks_file)
    
    output:
    tuple val(meta), path("*.chunk_*.vcf.gz"), path("*.chunk_*.vcf.gz.tbi"), emit: chunks
    path "versions.yml"                                                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Read chunks file and split VCF accordingly
    while IFS=\$'\\t' read -r chrom start end chunk_id; do
        echo "Processing chunk: \${chunk_id} (\${chrom}:\${start}-\${end})"
        
        # Extract chunk from VCF
        bcftools view \\
            --regions \${chrom}:\${start}-\${end} \\
            --output-type z \\
            --output ${prefix}.\${chunk_id}.vcf.gz \\
            ${vcf}
        
        # Index the chunk
        bcftools index -t ${prefix}.\${chunk_id}.vcf.gz
    done < ${chunks_file}
    
    # Count chunks created
    echo "Total chunks created: \$(ls ${prefix}.chunk_*.vcf.gz | wc -l)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chunk_0001.vcf.gz
    touch ${prefix}.chunk_0001.vcf.gz.tbi
    touch ${prefix}.chunk_0002.vcf.gz
    touch ${prefix}.chunk_0002.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}