process COMBINE_FREQ_BY_CHR {
    tag "${meta.id}_${chrm}"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(chrm), path(freq_files)
    
    output:
    tuple val(meta), val(chrm), path("*.chr_freq.tsv"), emit: chr_frequencies
    path "versions.yml"                                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    combine_freq_by_chr.py \\
        --prefix ${prefix} \\
        --chromosome ${chrm}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chrm}.chr_freq.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 1.5.0
        numpy: 1.23.0
    END_VERSIONS
    """
}