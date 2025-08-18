process FILTER_INFO_BY_TARGET {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_file)
    
    output:
    tuple val(meta), val(ref_name), path("*.filtered.info"), path("*.acc.info"), emit: filtered
    path "versions.yml"                                                         , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2_threshold = params.r2_threshold ?: 0.3
    """
    filter_info_by_target.py \\
        --input-file ${info_file} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --r2-threshold ${r2_threshold} \\
        ${args}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.filtered.info
    touch ${prefix}_${ref_name}.acc.info
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}