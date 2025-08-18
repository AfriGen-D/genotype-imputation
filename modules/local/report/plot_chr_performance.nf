process PLOT_CHR_PERFORMANCE {
    tag "${meta.dataset}_${meta.chromosome}"
    label 'python_plotting'
    publishDir "${params.outdir}/reports/chromosome/${meta.dataset}", mode: 'copy'
    
    input:
    tuple val(meta), path(chr_summary)
    
    output:
    tuple val(meta), path("*.chr_performance.pdf"), emit: plots
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.dataset}_${meta.chromosome}"
    """
    # Generate chromosome performance plots
    plot_chr_performance.py \\
        --chr-summary ${chr_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${meta.ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}