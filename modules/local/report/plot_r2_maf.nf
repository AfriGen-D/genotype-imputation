process PLOT_R2_MAF {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(well_info), path(acc_info)
    
    output:
    tuple val(meta), val(ref_name), path("*.r2_maf.png"), emit: plot
    path "versions.yml"                                 , emit: versions
    
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
    import numpy as np
    import sys
    
    # Read accuracy info file
    mafs = []
    rsqs = []
    
    with open("${acc_info}", 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\\t')
                if len(parts) >= 7:
                    try:
                        maf = float(parts[4])
                        rsq = float(parts[6])
                        mafs.append(maf)
                        rsqs.append(rsq)
                    except (ValueError, IndexError):
                        continue
    
    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if mafs and rsqs:
        # Create hexbin plot for large datasets
        if len(mafs) > 1000:
            hb = ax.hexbin(mafs, rsqs, gridsize=50, cmap='YlOrRd', mincnt=1)
            cb = plt.colorbar(hb)
            cb.set_label('Count')
        else:
            ax.scatter(mafs, rsqs, alpha=0.5, s=10)
        
        # Add reference lines
        ax.axhline(y=0.3, color='r', linestyle='--', alpha=0.5, label='R² = 0.3')
        ax.axhline(y=0.8, color='g', linestyle='--', alpha=0.5, label='R² = 0.8')
        
        ax.set_xlabel('Minor Allele Frequency (MAF)')
        ax.set_ylabel('Imputation Quality (R²)')
        ax.set_title(f'R² vs MAF - ${meta.id} (${ref_name})')
        ax.set_xlim([0, 0.5])
        ax.set_ylim([0, 1])
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig("${prefix}_${ref_name}.r2_maf.png", dpi=150)
    plt.close()
    
    print(f"R² vs MAF plot saved: ${prefix}_${ref_name}.r2_maf.png")
    print(f"Total variants plotted: {len(mafs)}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.r2_maf.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}