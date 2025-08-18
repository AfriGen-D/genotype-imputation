process COMBINE_CHUNK_REPORTS {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports_dataset/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(chunk_reports)
    
    output:
    tuple val(meta), val(ref_name), path("*.summary.json"), path("*.summary.txt"), emit: combined
    path "versions.yml"                                                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Use efficient streaming aggregation for large numbers of chunks
    aggregate_chunk_metrics.py \\
        --chunk-reports *.txt *.tsv *.json \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset-id ${meta.id} \\
        --max-chunks-detail 100 \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.summary.json
    touch ${prefix}_${ref_name}.summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}