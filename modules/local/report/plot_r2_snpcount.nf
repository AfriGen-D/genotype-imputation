process PLOT_R2_SNPCOUNT {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*.r2_snpcount.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plot_out = "${prefix}_${ref_name}_r2_SNPcount.pdf"
    def impute_info_cutoff = params.impute_info_cutoff ?: 0.3
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import glob
    
    # Combine all info files
    all_data = []
    for info_file in glob.glob("*.info*"):
        try:
            df = pd.read_csv(info_file, sep='\\t')
            if 'Rsq' in df.columns or 'R2' in df.columns:
                all_data.append(df)
        except:
            continue
    
    if not all_data:
        # Create empty plot if no data
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax.set_title('R² vs SNP Count')
        plt.savefig("${plot_out}")
    else:
        combined = pd.concat(all_data, ignore_index=True)
        
        # Get R² values
        if 'Rsq' in combined.columns:
            combined['r2'] = combined['Rsq']
        elif 'R2' in combined.columns:
            combined['r2'] = combined['R2']
        
        # Create R² bins
        r2_bins = np.arange(0, 1.05, 0.05)
        combined['r2_bin'] = pd.cut(combined['r2'], bins=r2_bins)
        
        # Count SNPs in each bin
        counts = combined.groupby('r2_bin').size().reset_index(name='count')
        
        # Calculate mean R² for each bin
        mean_r2 = combined.groupby('r2_bin')['r2'].mean().reset_index(name='mean_r2')
        
        # Merge counts and mean R²
        result = pd.merge(counts, mean_r2, on='r2_bin')
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Bar plot
        x_pos = np.arange(len(result))
        bars = ax.bar(x_pos, result['count'], alpha=0.7)
        
        # Color bars by R² threshold
        for i, (bar, r2) in enumerate(zip(bars, result['mean_r2'])):
            if r2 >= ${impute_info_cutoff}:
                bar.set_color('green')
            else:
                bar.set_color('red')
        
        # Labels
        ax.set_xlabel('R² Bin')
        ax.set_ylabel('Number of SNPs')
        ax.set_title(f'SNP Count by R² Bin\\n{meta.id} - {ref_name}')
        
        # X-axis labels
        labels = [f"{b.left:.2f}-{b.right:.2f}" for b in result['r2_bin']]
        ax.set_xticks(x_pos[::2])  # Show every other label to avoid crowding
        ax.set_xticklabels(labels[::2], rotation=45, ha='right')
        
        # Add threshold line
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax.text(0.02, 0.98, f"R² threshold: {impute_info_cutoff}", 
                transform=ax.transAxes, va='top', fontsize=10)
        
        plt.tight_layout()
        plt.savefig("${plot_out}")
    
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
    touch ${prefix}_${ref_name}_r2_SNPcount.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}