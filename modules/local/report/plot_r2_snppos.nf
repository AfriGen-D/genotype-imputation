process PLOT_R2_SNPPOS {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_file)
    
    output:
    tuple val(meta), val(ref_name), path("*.r2_snppos.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_${ref_name}_r2_SNPpos.pdf"
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Read info file
    df = pd.read_csv("${info_file}", sep='\\t')
    
    # Extract chromosome and position
    if 'SNP' in df.columns:
        df[['chr', 'pos']] = df['SNP'].str.split(':', expand=True)
        df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    elif 'CHROM' in df.columns and 'POS' in df.columns:
        df['chr'] = df['CHROM']
        df['pos'] = df['POS']
    
    # Get R² values
    if 'Rsq' in df.columns:
        df['r2'] = df['Rsq']
    elif 'R2' in df.columns:
        df['r2'] = df['R2']
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot R² vs position
    scatter = ax.scatter(df['pos'], df['r2'], alpha=0.5, s=1)
    
    # Add rolling mean
    window_size = max(1, len(df) // 100)
    df_sorted = df.sort_values('pos')
    rolling_mean = df_sorted['r2'].rolling(window=window_size, center=True).mean()
    ax.plot(df_sorted['pos'], rolling_mean, 'r-', linewidth=2, label='Rolling mean')
    
    ax.set_xlabel('Position')
    ax.set_ylabel('R²')
    ax.set_title(f'Imputation Quality (R²) vs SNP Position\\n{meta.id} - {ref_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("${output}")
    
    # Version info
    import sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
        f.write(f'    matplotlib: {plt.matplotlib.__version__}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}_r2_SNPpos.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}