process PLOT_CONCORDANCE_MAF {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(imputed_vcf), path(true_vcf), val(chr)
    
    output:
    tuple val(meta), val(ref_name), path("*_concordance_maf.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_suffix = chr ? "_chr${chr}" : "_genome"
    def plot_out = "${prefix}_${ref_name}${chr_suffix}_concordance_maf.pdf"
    def chr_arg = chr ? "--chr ${chr}" : ""
    def true_vcf_arg = true_vcf.name != 'NO_FILE' ? "--true-vcf ${true_vcf}" : ""
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/plot_concordance_maf.py .
    
    # Run the script with appropriate arguments
    python3 plot_concordance_maf.py \\
        ${imputed_vcf} \\
        ${plot_out} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}" \\
        ${chr_arg} \\
        ${true_vcf_arg}
    
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
    
    # Ensure file has proper permissions
    chmod 644 ${plot_out}
    
    # Small delay to ensure filesystem operations complete
    sleep 1
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_suffix = chr ? "_chr${chr}" : "_genome"
    """
    touch ${prefix}_${ref_name}${chr_suffix}_concordance_maf.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
        numpy: 1.24.0
    END_VERSIONS
    """
}