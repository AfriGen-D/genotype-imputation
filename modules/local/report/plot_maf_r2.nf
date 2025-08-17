process PLOT_MAF_R2 {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*_maf_r2.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plot_out = "${prefix}_${ref_name}_maf_r2.pdf"
    def impute_info_cutoff = params.impute_info_cutoff ?: 0.3
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/plot_maf_r2.py .
    
    # Run the script with appropriate arguments
    python3 plot_maf_r2.py \\
        ${plot_out} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}" \\
        --info-cutoff ${impute_info_cutoff}
    
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
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}_maf_r2.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}