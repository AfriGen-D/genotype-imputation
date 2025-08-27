process PLOT_AGGREGATED_R2_DASHBOARD {
    tag "$meta.id"
    label 'process_high'
    label 'python_plotting'
    
    container 'mamana/python-plotting:1.1.0'
    
    publishDir "${params.outdir}/reports/aggregated_quality/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), path("*_aggregated_dashboard.pdf")    , emit: plot
    tuple val(meta), path("*_aggregated_dashboard.png")    , emit: plot_png  
    tuple val(meta), path("*_aggregated_stats.json")       , emit: stats
    path "versions.yml"                                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def title = meta.id ? "${meta.id} - ${ref_name}" : ref_name
    """
    # Handle both .sites.vcf.gz and .info files
    if [ \$(ls *.sites.vcf.gz 2>/dev/null | wc -l) -gt 0 ]; then
        # VCF sites files
        python ${projectDir}/bin/plot_aggregated_r2_scores.py \\
            ${info_files} \\
            ${prefix} \\
            --title "${title}"
    elif [ \$(ls *.info 2>/dev/null | wc -l) -gt 0 ]; then
        # Legacy .info files  
        python ${projectDir}/bin/plot_aggregated_r2_scores.py \\
            ${info_files} \\
            ${prefix} \\
            --title "${title}"
    else
        echo "No info or sites.vcf.gz files found"
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aggregated_dashboard.pdf
    touch ${prefix}.aggregated_dashboard.png
    echo '{"summary": {"total_variants": 50000, "mean_r2": 0.82}}' > ${prefix}.aggregated_stats.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        matplotlib: 3.5.0
        pandas: 1.3.0
        numpy: 1.21.0
    END_VERSIONS
    """
}