process AVERAGE_R2 {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/rsquared/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*.average_r2.txt"), emit: average
    tuple val(meta), val(ref_name), path("*.r2_summary.csv"), emit: summary
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def meanr2_out = "${prefix}_${ref_name}.average_r2.txt"
    def summary_out = "${prefix}_${ref_name}.r2_summary.csv"
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/average_r2.py .
    
    # Run the script with appropriate arguments
    python3 average_r2.py \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}" \\
        --output-txt ${meanr2_out} \\
        --output-csv ${summary_out}
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Average R²: 0.85" > ${prefix}_${ref_name}.average_r2.txt
    echo "file,mean_r2,n_variants" > ${prefix}_${ref_name}.r2_summary.csv
    echo "test.info,0.85,1000" >> ${prefix}_${ref_name}.r2_summary.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        numpy: 1.24.0
    END_VERSIONS
    """
}