process PLOT_HIST_R2_SNPCOUNT {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*.r2_snpcount_hist.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plot_out = "${prefix}_${ref_name}_r2_SNPcount_hist.pdf"
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
        ax.set_title('R² Distribution Histogram')
        plt.savefig("${plot_out}")
    else:
        combined = pd.concat(all_data, ignore_index=True)
        
        # Get R² values
        if 'Rsq' in combined.columns:
            r2_values = combined['Rsq']
        elif 'R2' in combined.columns:
            r2_values = combined['R2']
        else:
            r2_values = pd.Series([])
        
        # Create histogram plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
        
        # Top plot: Full R² distribution
        n, bins, patches = ax1.hist(r2_values, bins=50, alpha=0.7, edgecolor='black')
        
        # Color bars by threshold
        for i, patch in enumerate(patches):
            if bins[i] >= ${impute_info_cutoff}:
                patch.set_facecolor('green')
            else:
                patch.set_facecolor('red')
        
        ax1.axvline(x=${impute_info_cutoff}, color='blue', linestyle='--', 
                    linewidth=2, label=f'Threshold: {impute_info_cutoff}')
        ax1.set_xlabel('R²')
        ax1.set_ylabel('Frequency')
        ax1.set_title(f'R² Distribution - All SNPs\\n{meta.id} - {ref_name}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Bottom plot: Cumulative distribution
        sorted_r2 = np.sort(r2_values)
        cumulative = np.arange(1, len(sorted_r2) + 1) / len(sorted_r2)
        
        ax2.plot(sorted_r2, cumulative, linewidth=2)
        ax2.axvline(x=${impute_info_cutoff}, color='blue', linestyle='--', 
                    linewidth=2, label=f'Threshold: {impute_info_cutoff}')
        ax2.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
        ax2.set_xlabel('R²')
        ax2.set_ylabel('Cumulative Proportion')
        ax2.set_title('Cumulative R² Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Add statistics text
        well_imputed = (r2_values >= ${impute_info_cutoff}).sum()
        total = len(r2_values)
        pct_well = (well_imputed / total * 100) if total > 0 else 0
        
        stats_text = f"Total SNPs: {total:,}\\n"
        stats_text += f"Well imputed (R² ≥ {impute_info_cutoff}): {well_imputed:,} ({pct_well:.1f}%)\\n"
        stats_text += f"Mean R²: {r2_values.mean():.3f}\\n"
        stats_text += f"Median R²: {r2_values.median():.3f}"
        
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
                 va='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
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
    touch ${prefix}_${ref_name}_r2_SNPcount_hist.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}