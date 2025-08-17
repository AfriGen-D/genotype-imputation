process PLOT_R2_SNPPOS {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_file)
    
    output:
    tuple val(meta), val(ref_name), path("*_r2_snppos.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_${ref_name}_r2_snppos.pdf"
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/plot_r2_snppos.py .
    
    # Run the script with appropriate arguments
    python3 plot_r2_snppos.py \\
        ${info_file} \\
        ${output} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}"
    
    # Ensure filesystem sync
    sync
    
    # Verify file exists
    if [ ! -f "${output}" ]; then
        echo "ERROR: Expected output file ${output} was not created!"
        ls -la *.pdf || echo "No PDF files found"
        exit 1
    fi
    
    echo "Output file verified: ${output}"
    ls -la ${output}
    
    # Ensure file has proper permissions
    chmod 644 ${output}
    
    # Double-check file is in current directory
    echo "Current directory contents:"
    pwd
    ls -la *.pdf 2>/dev/null || true
    
    # Small delay to ensure filesystem operations complete
    sleep 1
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}_r2_snppos.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}