process EXTRACT_IMPUTE_INFO {
    tag "$meta.id"
    label 'process_single'
    
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"
    
    input:
    tuple val(meta), val(ref_name), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), val(ref_name), path("*.info.gz"), emit: info
    path "versions.yml"                              , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import gzip
    import sys
    
    # Extract info from Minimac4 output
    info_file = "${prefix}_${ref_name}.info.gz"
    
    # For now, create a simple info file
    # In production, this would parse the actual Minimac4 info output
    with gzip.open(info_file, 'wt') as f:
        f.write("SNP\\tREF(0)\\tALT(1)\\tALT_Frq\\tMAF\\tAvgCall\\tRsq\\tGenotyped\\tLooRsq\\tEmpR\\tEmpRsq\\tDose0\\tDose1\\n")
        # Add placeholder data
        f.write("chr1:1000\\tA\\tG\\t0.15\\t0.15\\t0.99\\t0.95\\tGenotyped\\t0.94\\t0.96\\t0.95\\t0.01\\t0.02\\n")
    
    print(f"Info file created: {info_file}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.info.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}