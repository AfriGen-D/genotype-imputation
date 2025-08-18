process SUMMARIZE_DATASET_METRICS {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports_dataset/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(summary_json), path(summary_txt)
    
    output:
    tuple val(meta), val(ref_name), path("*.dataset_summary.json"), emit: summary
    path "versions.yml"                                           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # The aggregation already created the summary, just rename for clarity
    cp ${summary_json} ${prefix}_${ref_name}.dataset_summary.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.dataset_summary.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.0.0
    END_VERSIONS
    """
}