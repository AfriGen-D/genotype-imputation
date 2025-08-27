process GENERATE_SUMMARY_REPORT {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(well_summary), path(accuracy_txt), path(avg_r2)
    
    output:
    tuple val(meta), val(ref_name), path("*_summary_report.html"), emit: report
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def report_out = "${prefix}_${ref_name}_summary_report.html"
    def metrics_dir = metrics_files.name != 'NO_FILES' ? "--metrics-dir ." : ""
    def plots_dir = plot_files.name != 'NO_FILES' ? "--plots-dir ." : ""
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/generate_summary_report.py .
    
    # Run the script with appropriate arguments
    python3 generate_summary_report.py \\
        ${report_out} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}" \\
        ${metrics_dir} \\
        ${plots_dir}
    
    # Ensure filesystem sync
    sync
    
    # Verify file exists
    if [ ! -f "${report_out}" ]; then
        echo "ERROR: Expected output file ${report_out} was not created!"
        ls -la *.html || echo "No HTML files found"
        exit 1
    fi
    
    echo "Output file verified: ${report_out}"
    ls -la ${report_out}
    
    # Ensure file has proper permissions
    chmod 644 ${report_out}
    
    # Small delay to ensure filesystem operations complete
    sleep 1
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<'EOF' > ${prefix}_${ref_name}_summary_report.html
    <!DOCTYPE html>
    <html>
    <head><title>Summary Report</title></head>
    <body><h1>Imputation Summary Report</h1><p>Stub report for ${meta.id}</p></body>
    </html>
    EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        numpy: 1.24.0
    END_VERSIONS
    """
}