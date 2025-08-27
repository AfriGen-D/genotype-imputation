process PLOT_IMPUTATION_ACCURACY_MAF_BINS {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'
    
    container 'mamana/python-plotting:1.1.0'
    
    publishDir "${params.outdir}/reports/imputation_quality/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_file)
    
    output:
    tuple val(meta), path("*_maf_accuracy.pdf")        , emit: plot
    tuple val(meta), path("*_maf_accuracy.png")        , emit: plot_png
    tuple val(meta), path("*_stats.json")              , emit: stats
    tuple val(meta), path("*_summary.txt")             , emit: summary
    path "versions.yml"                                 , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = meta.id ?: "Sample"
    """
    # Run MAF bins accuracy analysis
    python ${projectDir}/bin/plot_imputation_accuracy_maf_bins.py \\
        ${info_file} \\
        ${prefix} \\
        --title "Imputation Accuracy Analysis - ${ref_name}" \\
        --sample-name "${sample_name}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf_accuracy.pdf
    touch ${prefix}.maf_accuracy.png
    echo '{"total_variants": 1000, "mean_r2": 0.85}' > ${prefix}.stats.json
    echo "Summary statistics" > ${prefix}.summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        matplotlib: 3.5.0
        pandas: 1.3.0
    END_VERSIONS
    """
}