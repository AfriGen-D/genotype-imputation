process COMBINE_FREQ_GENOME {
    tag "$meta.id"
    label 'process_high'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), path(chr_freq_files)
    
    output:
    tuple val(meta), path("*.genome_freq.tsv")     , emit: genome_frequencies
    tuple val(meta), path("*.genome_freq_summary.txt"), emit: summary
    path "versions.yml"                            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    combine_freq_genome.py \\
        --output-prefix ${prefix} \\
        --input-pattern "*.chr_freq.tsv" \\
        ${args}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genome_freq.tsv
    touch ${prefix}.genome_freq_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 1.5.0
        numpy: 1.23.0
    END_VERSIONS
    """
}