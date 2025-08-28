process GENERATE_REPORT_INDEX {
    tag "$dataset"
    label 'process_single'
    container 'mamana/python-plotting:1.1.0'
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true
    
    input:
    val(dataset)
    path(reports_dir)
    
    output:
    path("${dataset}_report_index.html"), emit: index
    path("report_index.html"), emit: main_index
    path "versions.yml", emit: versions
    
    script:
    """
    # Create main index for all reports
    python3 ${projectDir}/bin/generate_report_index.py \
        ${reports_dir} \
        --output report_index.html \
        --dataset "${dataset}"
    
    # Create dataset-specific index
    python3 ${projectDir}/bin/generate_report_index.py \
        ${reports_dir} \
        --output "${dataset}_report_index.html" \
        --dataset "${dataset}"
    
    # Create version file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
    
    stub:
    """
    touch "${dataset}_report_index.html"
    touch "report_index.html"
    touch versions.yml
    """
}