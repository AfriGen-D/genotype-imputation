process PLOT_DATASET_SUMMARY {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots_dataset/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(summary_json), path(chunk_plots)
    
    output:
    tuple val(meta), val(ref_name), path("*.dataset_summary.png"), emit: plot
    path "versions.yml"                                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_dataset_summary.py \\
        --summary-json ${summary_json} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.dataset_summary.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}