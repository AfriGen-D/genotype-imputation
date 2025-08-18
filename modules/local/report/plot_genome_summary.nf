process PLOT_GENOME_SUMMARY {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/genome/${meta.dataset}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(genome_summary), path(chr_plots)
    
    output:
    tuple val(meta), val(ref_name), path("*.genome_plots.pdf"), emit: genome_plots
    path "versions.yml"                                        , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_plots_arg = chr_plots ? "--chr-plots ${chr_plots}" : ""
    """
    # Generate genome-wide plots
    plot_genome_summary.py \\
        --genome-summary ${genome_summary} \\
        ${chr_plots_arg} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.genome_plots.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}