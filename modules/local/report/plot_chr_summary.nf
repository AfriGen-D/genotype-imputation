process PLOT_CHR_SUMMARY {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/chromosome/${meta.dataset}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(chr_summary), path(chunk_plots)
    
    output:
    tuple val(meta), val(ref_name), path("*.chr_plots.pdf"), emit: chr_plots
    path "versions.yml"                                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_plots_arg = chunk_plots ? "--chunk-plots ${chunk_plots}" : ""
    """
    # Generate chromosome-level plots
    plot_chromosome_summary.py \\
        --chr-summary ${chr_summary} \\
        ${chunk_plots_arg} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.chr_plots.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}