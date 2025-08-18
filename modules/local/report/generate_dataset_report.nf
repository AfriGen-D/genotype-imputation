process GENERATE_DATASET_REPORT {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/final/${meta.dataset}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(genome_summary), path(genome_plots)
    
    output:
    tuple val(meta), val(ref_name), path("*.final_report.html"), emit: report
    tuple val(meta), val(ref_name), path("*.final_report.pdf") , emit: report_pdf
    path "versions.yml"                                         , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Generate comprehensive HTML and PDF report
    generate_final_report.py \\
        --genome-summary ${genome_summary} \\
        --genome-plots ${genome_plots} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --population ${meta.population} \\
        --study ${meta.study} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.final_report.html
    touch ${prefix}_${ref_name}.final_report.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}