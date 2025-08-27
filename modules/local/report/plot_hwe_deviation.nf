process PLOT_HWE_DEVIATION {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(vcf_file), path(vcf_index)
    
    output:
    tuple val(meta), val(ref_name), path("*_hwe_deviation.pdf"), emit: plot
    tuple val(meta), val(ref_name), path("*_hwe_deviation_significant.txt"), optional: true, emit: significant
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plot_out = "${prefix}_${ref_name}_hwe_deviation.pdf"
    def p_threshold = params.hwe_p_threshold ?: 1e-6
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/plot_hwe_deviation.py .
    
    # Run the script with appropriate arguments
    python3 plot_hwe_deviation.py \\
        ${vcf_file} \\
        ${plot_out} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}" \\
        --p-threshold ${p_threshold}
    
    # Ensure filesystem sync
    sync
    
    # Verify file exists
    if [ ! -f "${plot_out}" ]; then
        echo "ERROR: Expected output file ${plot_out} was not created!"
        ls -la *.pdf || echo "No PDF files found"
        exit 1
    fi
    
    echo "Output file verified: ${plot_out}"
    ls -la ${plot_out}
    
    # Check for significant variants file
    sig_file="${plot_out%.pdf}_significant.txt"
    if [ -f "\${sig_file}" ]; then
        echo "Significant variants file found: \${sig_file}"
        ls -la "\${sig_file}"
    fi
    
    # Ensure file has proper permissions
    chmod 644 ${plot_out}
    [ -f "\${sig_file}" ] && chmod 644 "\${sig_file}"
    
    # Small delay to ensure filesystem operations complete
    sleep 1
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_suffix = chr ? "_chr${chr}" : "_genome"
    """
    touch ${prefix}_${ref_name}${chr_suffix}_hwe_deviation.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
        numpy: 1.24.0
        scipy: 1.5.0
    END_VERSIONS
    """
}