process CHECK_OVERLAP {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), path(chunk_map), path(ref_map)
    
    output:
    tuple val(meta), path("*.overlap.txt"), path("*.overlap.status"), emit: overlap
    path "versions.yml"                                              , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minRatio = params.minRatio ?: 0.001
    def chrm = meta.contig ?: 'chr21'
    def chunk_start = meta.start ?: ''
    def chunk_end = meta.end ?: ''
    def region = "${chrm}:${chunk_start}-${chunk_end}"
    """
    # Calculate overlap using Python script
    python3 /users/mamana/genotype-imputation/bin/calculate_overlap.py \\
        --chunk-map ${chunk_map} \\
        --ref-map ${ref_map} \\
        --min-ratio ${minRatio} \\
        --output ${prefix}.overlap.txt \\
        --chunk-id ${prefix} \\
        --region "${region}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.overlap.txt
    echo "PASS" > ${prefix}.overlap.status
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}