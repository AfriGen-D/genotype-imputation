process PLOT_CHR_ALL_METRICS {
    tag "${meta.dataset}_${meta.chromosome}"
    label 'python_plotting'
    publishDir "${params.outdir}/reports/chromosome/${meta.dataset}/plots", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(chr_summary), path(chunk_plots)
    
    output:
    tuple val(meta), path("*.chr_performance.pdf"), emit: performance
    tuple val(meta), path("*.chr_r2_position.pdf"), emit: r2_position
    tuple val(meta), path("*.chr_accuracy.pdf"), emit: accuracy
    tuple val(meta), path("*.chr_maf_analysis.pdf"), emit: maf_analysis
    tuple val(meta), path("*.chr_freq_comparison.pdf"), emit: freq_comparison
    tuple val(meta), path("*.chr_all_metrics.pdf"), emit: combined
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.dataset}_${meta.chromosome}"
    """
    # Generate all chromosome-level plots
    
    # 1. Performance metrics
    plot_chr_performance.py \\
        --chr-summary ${chr_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    # 2. R² by position
    plot_chr_r2_position.py \\
        --chr-summary ${chr_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    # 3. Accuracy metrics
    plot_chr_accuracy.py \\
        --chr-summary ${chr_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    # 4. MAF analysis
    plot_chr_maf_analysis.py \\
        --chr-summary ${chr_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    # 5. Frequency comparison
    plot_chr_freq_comparison.py \\
        --chr-summary ${chr_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    # 6. Combined comprehensive report
    combine_chr_plots.py \\
        --chr-summary ${chr_summary} \\
        --performance-plot ${prefix}_${ref_name}.chr_performance.pdf \\
        --r2-position-plot ${prefix}_${ref_name}.chr_r2_position.pdf \\
        --accuracy-plot ${prefix}_${ref_name}.chr_accuracy.pdf \\
        --maf-plot ${prefix}_${ref_name}.chr_maf_analysis.pdf \\
        --freq-plot ${prefix}_${ref_name}.chr_freq_comparison.pdf \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}