process COMBINE_INFO {
    tag "$meta.id"
    label 'process_single'
    
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*.combined.info.gz"), emit: info
    path "versions.yml"                                        , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import gzip
    import sys
    from pathlib import Path
    
    info_files = "${info_files}".split()
    output_file = "${prefix}_${ref_name}.combined.info.gz"
    
    # Combine all info files
    with gzip.open(output_file, 'wt') as outf:
        header_written = False
        
        for i, info_file in enumerate(sorted(info_files)):
            with gzip.open(info_file, 'rt') if info_file.endswith('.gz') else open(info_file, 'r') as inf:
                lines = inf.readlines()
                
                # Write header only once
                if not header_written and lines and lines[0].startswith('SNP'):
                    outf.write(lines[0])
                    header_written = True
                    start_idx = 1
                else:
                    start_idx = 1 if lines and lines[0].startswith('SNP') else 0
                
                # Write data lines
                for line in lines[start_idx:]:
                    if line.strip():
                        outf.write(line)
    
    print(f"Combined {len(info_files)} info files into {output_file}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.combined.info.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}