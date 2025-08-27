process COMPARE_PRE_POST_IMPUTATION {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.1.0'
    
    publishDir "${params.outdir}/reports/pre_post_comparison/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), path(pre_vcf), path(pre_index), path(post_vcf), path(post_index)
    
    output:
    tuple val(meta), path("*.comparison_stats.txt")     , emit: stats
    tuple val(meta), path("*.comparison_stats.txt")     , emit: comparison  // Alias for compatibility
    tuple val(meta), path("*.variant_counts.png")       , emit: count_plot
    tuple val(meta), path("*.maf_distribution.png")     , emit: maf_plot
    tuple val(meta), path("*.coverage_improvement.png") , emit: coverage_plot
    tuple val(meta), path("*.variant_gain.png")         , emit: gain_plot
    tuple val(meta), path("*.comparison_report.html")   , emit: report
    path "versions.yml"                                  , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # For now, create placeholder outputs to test module compilation
    echo "Pre-imputation VCF: ${pre_vcf}" > ${prefix}.comparison_stats.txt
    echo "Post-imputation VCF: ${post_vcf}" >> ${prefix}.comparison_stats.txt
    
    # Create placeholder files
    touch ${prefix}.variant_counts.png
    touch ${prefix}.maf_distribution.png
    touch ${prefix}.coverage_improvement.png
    touch ${prefix}.variant_gain.png
    
    # Create simple HTML report
    cat <<EOF > ${prefix}.comparison_report.html
    <!DOCTYPE html>
    <html>
    <head><title>Comparison Report</title></head>
    <body>
    <h1>Pre/Post Imputation Comparison</h1>
    <p>Sample: ${prefix}</p>
    <p>Pre-imputation: ${pre_vcf}</p>
    <p>Post-imputation: ${post_vcf}</p>
    </body>
    </html>
    EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}