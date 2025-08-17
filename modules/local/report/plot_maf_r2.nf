process PLOT_MAF_R2 {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/plots/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*.maf_r2.pdf"), emit: plot
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plot_out = "${prefix}_${ref_name}_MAF_r2.pdf"
    def impute_info_cutoff = params.impute_info_cutoff ?: 0.3
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/plot_maf_r2.py .
    
    # Run the script with appropriate arguments
    python3 plot_maf_r2.py \\
        ${plot_out} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}" \\
        --info-cutoff ${impute_info_cutoff}
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import glob
    
    # Combine all info files
    all_data = []
    for info_file in glob.glob("*.info*"):
        try:
            df = pd.read_csv(info_file, sep='\\t')
            if ('Rsq' in df.columns or 'R2' in df.columns) and ('MAF' in df.columns or 'AF' in df.columns):
                all_data.append(df)
        except:
            continue
    
    if not all_data:
        # Create empty plot if no data
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax.set_title('MAF vs R²')
        plt.savefig("${plot_out}")
    else:
        combined = pd.concat(all_data, ignore_index=True)
        
        # Get R² values
        if 'Rsq' in combined.columns:
            combined['r2'] = combined['Rsq']
        elif 'R2' in combined.columns:
            combined['r2'] = combined['R2']
        
        # Get MAF values
        if 'MAF' in combined.columns:
            combined['maf'] = combined['MAF']
        elif 'AF' in combined.columns:
            # Convert AF to MAF
            combined['maf'] = combined['AF'].apply(lambda x: min(x, 1-x) if pd.notna(x) else x)
        
        # Remove missing values
        combined = combined.dropna(subset=['r2', 'maf'])
        
        # Create MAF bins
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        combined['maf_bin'] = pd.cut(combined['maf'], bins=maf_bins)
        
        # Calculate statistics for each MAF bin
        stats = combined.groupby('maf_bin').agg({
            'r2': ['mean', 'median', 'std', 'count']
        }).round(3)
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Top plot: Scatter plot of MAF vs R²
        scatter = ax1.scatter(combined['maf'], combined['r2'], 
                             alpha=0.3, s=1, c=combined['r2'], 
                             cmap='RdYlGn', vmin=0, vmax=1)
        
        # Add mean R² line for each MAF bin
        for maf_bin in stats.index:
            if pd.notna(maf_bin):
                bin_data = combined[combined['maf_bin'] == maf_bin]
                if len(bin_data) > 0:
                    maf_center = (maf_bin.left + maf_bin.right) / 2
                    mean_r2 = stats.loc[maf_bin, ('r2', 'mean')]
                    ax1.plot(maf_center, mean_r2, 'ko', markersize=8)
        
        ax1.axhline(y=${impute_info_cutoff}, color='blue', linestyle='--', 
                    linewidth=1, label=f'R² threshold: {impute_info_cutoff}')
        ax1.set_xlabel('Minor Allele Frequency (MAF)')
        ax1.set_ylabel('R²')
        ax1.set_title(f'MAF vs Imputation Quality (R²)\\n{meta.id} - {ref_name}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax1, label='R²')
        
        # Bottom plot: Box plot by MAF bins
        maf_labels = [f"{b.left:.2f}-{b.right:.2f}" for b in stats.index if pd.notna(b)]
        box_data = [combined[combined['maf_bin'] == b]['r2'].values 
                    for b in stats.index if pd.notna(b)]
        
        bp = ax2.boxplot(box_data, labels=maf_labels, patch_artist=True)
        
        # Color boxes by mean R²
        for patch, maf_bin in zip(bp['boxes'], stats.index):
            if pd.notna(maf_bin):
                mean_r2 = stats.loc[maf_bin, ('r2', 'mean')]
                if mean_r2 >= ${impute_info_cutoff}:
                    patch.set_facecolor('lightgreen')
                else:
                    patch.set_facecolor('lightcoral')
        
        ax2.axhline(y=${impute_info_cutoff}, color='blue', linestyle='--', 
                    linewidth=1, label=f'R² threshold: {impute_info_cutoff}')
        ax2.set_xlabel('MAF Bins')
        ax2.set_ylabel('R²')
        ax2.set_title('R² Distribution by MAF Bins')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Add statistics table
        stats_text = "MAF Bin Statistics:\\n"
        for i, (maf_bin, row) in enumerate(stats.iterrows()):
            if pd.notna(maf_bin) and i < 5:  # Show first 5 bins
                stats_text += f"{maf_bin.left:.2f}-{maf_bin.right:.2f}: "
                stats_text += f"mean={row[('r2', 'mean')]:.3f}, n={int(row[('r2', 'count')])}\\n"
        
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
                 va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
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
    touch ${prefix}_${ref_name}_MAF_r2.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}