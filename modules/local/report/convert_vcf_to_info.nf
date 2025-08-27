process CONVERT_VCF_TO_INFO {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.1.0'
    
    input:
    tuple val(meta), val(ref_name), path(sites_vcf)
    
    output:
    tuple val(meta), val(ref_name), path("*.info"), emit: info
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = "${prefix}_${ref_name}.info"
    """
    python ${projectDir}/bin/vcf_to_info.py \\
        ${sites_vcf} \\
        ${output_file}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.info
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}