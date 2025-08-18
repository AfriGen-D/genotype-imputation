process PLOT_R2_MAF {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(well_info), path(acc_info)
    
    output:
    tuple val(meta), val(ref_name), path("*.r2_maf.png"), emit: plot
    path "versions.yml"                                 , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_r2_maf.py \\
        --acc-info-file ${acc_info} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --sample-id ${meta.id} \\
        ${args}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.r2_maf.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}