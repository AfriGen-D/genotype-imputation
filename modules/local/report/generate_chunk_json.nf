process GENERATE_CHUNK_JSON {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(accuracy_txt), path(well_imputed_txt), path(summary_txt)
    
    output:
    tuple val(meta), val(ref_name), path("*.chunk_summary.json"), emit: json
    path "versions.yml"                                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_chunk_json.py \\
        --accuracy-file ${accuracy_txt} \\
        --well-imputed-file ${well_imputed_txt} \\
        --summary-file ${summary_txt} \\
        --chunk-id ${meta.id} \\
        --ref-name ${ref_name} \\
        --output ${prefix}_${ref_name}.chunk_summary.json \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.chunk_summary.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}