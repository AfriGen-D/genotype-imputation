process PLOT_ACCURACY {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(accuracy_report), path(accuracy_tsv)
    
    output:
    tuple val(meta), val(ref_name), path("*.accuracy.png"), emit: plot
    path "versions.yml"                                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_accuracy.py \\
        --accuracy-tsv ${accuracy_tsv} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --sample-id ${meta.id} \\
        ${args}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.accuracy.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
        pandas: 2.0.0
    END_VERSIONS
    """
}