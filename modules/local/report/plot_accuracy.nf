process PLOT_ACCURACY {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(accuracy_report), path(accuracy_tsv)
    
    output:
    tuple val(meta), val(ref_name), path("*.accuracy.png"), emit: plot
    path "versions.yml"                                   , emit: versions
    
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
    
    # Read accuracy TSV
    data = pd.read_csv("${accuracy_tsv}", sep='\\t')
    
    # Create accuracy plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Mean Rsq by MAF
    if 'MAF_BIN' in data.columns and 'MEAN_RSQ' in data.columns:
        ax1.plot(range(len(data)), data['MEAN_RSQ'], 'o-', linewidth=2, markersize=8)
        ax1.set_xticks(range(len(data)))
        ax1.set_xticklabels(data['MAF_BIN'], rotation=45)
        ax1.set_xlabel('MAF Bin')
        ax1.set_ylabel('Mean Rsq')
        ax1.set_title(f'Imputation Accuracy by MAF')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([0, 1])
    
    # Plot 2: Count distribution
    if 'MAF_BIN' in data.columns and 'COUNT' in data.columns:
        ax2.bar(range(len(data)), data['COUNT'], alpha=0.7)
        ax2.set_xticks(range(len(data)))
        ax2.set_xticklabels(data['MAF_BIN'], rotation=45)
        ax2.set_xlabel('MAF Bin')
        ax2.set_ylabel('Number of Variants')
        ax2.set_title(f'Variant Distribution by MAF')
        ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f'Imputation Accuracy - ${meta.id} (${ref_name})', fontsize=14)
    plt.tight_layout()
    plt.savefig("${prefix}_${ref_name}.accuracy.png", dpi=150)
    plt.close()
    
    print(f"Accuracy plot saved: ${prefix}_${ref_name}.accuracy.png")
    
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
    touch ${prefix}_${ref_name}.accuracy.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
        pandas: 2.0.0
    END_VERSIONS
    """
}