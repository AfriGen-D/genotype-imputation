process PLOT_PERFORMANCE {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(well_imputed), path(summary)
    
    output:
    tuple val(meta), val(ref_name), path("*.performance.png"), emit: plot
    path "versions.yml"                                       , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    import pandas as pd
    import sys
    
    # Read the well imputed report
    data = pd.read_csv("${well_imputed}", sep='\\t')
    
    # Create performance plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if 'MAF_BIN' in data.columns and 'COUNT' in data.columns:
        ax.bar(range(len(data)), data['COUNT'])
        ax.set_xticks(range(len(data)))
        ax.set_xticklabels(data['MAF_BIN'], rotation=45)
        ax.set_xlabel('MAF Bin')
        ax.set_ylabel('Number of Well Imputed Variants')
        ax.set_title(f'Imputation Performance by MAF - ${meta.id}')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("${prefix}_${ref_name}.performance.png", dpi=150)
    plt.close()
    
    print(f"Performance plot saved: ${prefix}_${ref_name}.performance.png")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.performance.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
        pandas: 2.0.0
    END_VERSIONS
    """
}