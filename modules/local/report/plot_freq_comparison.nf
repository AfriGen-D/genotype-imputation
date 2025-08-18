process PLOT_FREQ_COMPARISON {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(level), path(freq_file)
    
    output:
    tuple val(meta), val(level), path("*.freq_comparison.png"), emit: plot
    path "*.freq_stats.txt"                                   , emit: stats
    path "versions.yml"                                       , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_freq_comparison.py \\
        ${freq_file} \\
        --prefix ${prefix} \\
        --level ${level}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${level}.freq_comparison.png
    touch ${prefix}_${level}.freq_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
        pandas: 1.5.0
        numpy: 1.23.0
        scipy: 1.9.0
    END_VERSIONS
    """
}