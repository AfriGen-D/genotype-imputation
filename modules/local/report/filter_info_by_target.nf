process FILTER_INFO_BY_TARGET {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(info_file)
    
    output:
    tuple val(meta), val(ref_name), path("*.filtered.info"), path("*.acc.info"), emit: filtered
    path "versions.yml"                                                         , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2_threshold = params.r2_threshold ?: 0.3
    """
    #!/usr/bin/env python3
    
    import gzip
    import sys
    
    input_file = "${info_file}"
    well_imputed_file = "${prefix}_${ref_name}.filtered.info"
    accuracy_file = "${prefix}_${ref_name}.acc.info"
    
    # Read and filter info file
    well_imputed = []
    accuracy = []
    
    with gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\\t')
                if len(parts) >= 7:
                    try:
                        rsq = float(parts[6])  # Assuming Rsq is in column 7
                        if rsq >= ${r2_threshold}:
                            well_imputed.append(line)
                        accuracy.append(line)
                    except (ValueError, IndexError):
                        continue
    
    # Write filtered files
    with open(well_imputed_file, 'w') as f:
        f.write(header)
        f.writelines(well_imputed)
    
    with open(accuracy_file, 'w') as f:
        f.write(header)
        f.writelines(accuracy)
    
    print(f"Well imputed variants (Rsq >= ${r2_threshold}): {len(well_imputed)}")
    print(f"Total variants for accuracy: {len(accuracy)}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.filtered.info
    touch ${prefix}_${ref_name}.acc.info
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}