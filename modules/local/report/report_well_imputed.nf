process REPORT_WELL_IMPUTED {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(well_info), path(acc_info)
    
    output:
    tuple val(meta), val(ref_name), path("*.well_imputed.txt"), path("*.well_imputed_summary.txt"), emit: report
    path "versions.yml"                                                                            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    report_well_imputed.py \\
        --well-info-file ${well_info} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --sample-id ${meta.id} \\
        ${args}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.well_imputed.txt
    touch ${prefix}_${ref_name}.well_imputed_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}