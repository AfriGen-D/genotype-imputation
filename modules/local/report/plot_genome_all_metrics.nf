process PLOT_GENOME_ALL_METRICS {
    tag "${meta.dataset}"
    label 'python_plotting'
    publishDir "${params.outdir}/reports/genome/${meta.dataset}/plots", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(genome_summary), path(chr_plots)
    
    output:
    tuple val(meta), val(ref_name), path("*.genome_performance.pdf"), emit: performance
    tuple val(meta), val(ref_name), path("*.genome_r2_distribution.pdf"), emit: r2_distribution
    tuple val(meta), val(ref_name), path("*.genome_accuracy.pdf"), emit: accuracy
    tuple val(meta), val(ref_name), path("*.genome_maf_analysis.pdf"), emit: maf_analysis
    tuple val(meta), val(ref_name), path("*.genome_freq_comparison.pdf"), emit: freq_comparison
    tuple val(meta), val(ref_name), path("*.genome_chr_comparison.pdf"), emit: chr_comparison
    tuple val(meta), val(ref_name), path("*.genome_all_metrics.pdf"), emit: combined
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.dataset}"
    """
    # Generate all genome-wide plots
    
    # 1. Overall performance metrics
    plot_genome_performance.py \\
        --genome-summary ${genome_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    # 2. R² distribution across genome
    plot_genome_r2_distribution.py \\
        --genome-summary ${genome_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    # 3. Accuracy metrics genome-wide
    plot_genome_accuracy.py \\
        --genome-summary ${genome_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    # 4. MAF analysis genome-wide
    plot_genome_maf_analysis.py \\
        --genome-summary ${genome_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    # 5. Frequency comparison genome-wide
    plot_genome_freq_comparison.py \\
        --genome-summary ${genome_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    # 6. Chromosome comparison
    plot_genome_chr_comparison.py \\
        --genome-summary ${genome_summary} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    # 7. Combined comprehensive report
    combine_genome_plots.py \\
        --genome-summary ${genome_summary} \\
        --performance-plot ${prefix}_${ref_name}.genome_performance.pdf \\
        --r2-plot ${prefix}_${ref_name}.genome_r2_distribution.pdf \\
        --accuracy-plot ${prefix}_${ref_name}.genome_accuracy.pdf \\
        --maf-plot ${prefix}_${ref_name}.genome_maf_analysis.pdf \\
        --freq-plot ${prefix}_${ref_name}.genome_freq_comparison.pdf \\
        --chr-plot ${prefix}_${ref_name}.genome_chr_comparison.pdf \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}