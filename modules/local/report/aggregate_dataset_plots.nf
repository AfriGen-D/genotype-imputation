process AGGREGATE_DATASET_PLOTS {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots_dataset/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(chunk_plots)
    
    output:
    tuple val(meta), val(ref_name), path("*.dataset_plots.pdf"), emit: aggregated_plots
    path "versions.yml"                                         , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    aggregate_dataset_plots.py \\
        --chunk-plots *.pdf *.png \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset-id ${meta.id} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.dataset_plots.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}