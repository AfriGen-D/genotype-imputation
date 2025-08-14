process PLOT_INFO_DISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'

    input:
    tuple val(meta), path(info_file)

    output:
    tuple val(meta), path("*.info_distribution.pdf"), emit: plot
    tuple val(meta), path("*.info_distribution_stats.tsv"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    import gzip
    from pathlib import Path
    from scipy import stats as scipy_stats

    # Set plotting style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")

    def read_file(filename):
        \"\"\"Read potentially gzipped files\"\"\"
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as f:
                return pd.read_csv(f, sep='\\t')
        else:
            return pd.read_csv(filename, sep='\\t')

    # Input parameters
    info_file = '${info_file}'
    output_prefix = '${prefix}'

    try:
        # Read the data
        data = read_file(info_file)
        print(f"Loaded {len(data)} variants from info file")
        print(f"Available columns: {list(data.columns)}")
        
        # Find INFO score column (various naming conventions)
        info_cols = ['INFO', 'AvgCall', 'Rsq', 'LooRsq', 'EmpRsq', 'R2']
        info_col = None
        
        for col in info_cols:
            if col in data.columns:
                info_col = col
                break
        
        if not info_col:
            raise ValueError(f"Could not find INFO score column. Available: {list(data.columns)}")
            
        print(f"Using INFO score column: {info_col}")
        
        # Clean the data
        clean_data = data[data[info_col] != '-'].copy()
        clean_data = clean_data[clean_data[info_col].notna()]
        
        # Convert to numeric
        clean_data['info_score'] = pd.to_numeric(clean_data[info_col], errors='coerce')
        clean_data = clean_data.dropna(subset=['info_score'])
        
        # Filter extreme values (INFO scores should be between 0 and 1)
        clean_data = clean_data[
            (clean_data['info_score'] >= 0) & 
            (clean_data['info_score'] <= 1)
        ]
        
        print(f"After cleaning: {len(clean_data)} variants with valid INFO scores")
        
        if len(clean_data) == 0:
            raise ValueError("No valid INFO scores found after filtering")
            
        info_scores = clean_data['info_score'].values
        
    except Exception as e:
        print(f"Error processing data: {e}")
        # Create error plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f'Error loading or processing data\\n{str(e)}', 
                ha='center', va='center', fontsize=12, transform=ax.transAxes)
        ax.set_title('INFO Score Distribution - Error', fontsize=14)
        plt.savefig(f'{output_prefix}.info_distribution.pdf', dpi=150, bbox_inches='tight')
        plt.close()
        
        # Create error stats
        with open(f'{output_prefix}.info_distribution_stats.tsv', 'w') as f:
            f.write('Statistic\\tValue\\n')
            f.write(f'Error\\t{str(e)}\\n')
        
        exit(1)

    # Create the visualization
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Main histogram (top row, spanning 2 columns)
    ax_main = fig.add_subplot(gs[0, :2])
    
    # Create histogram with multiple binning strategies
    bins = min(50, int(np.sqrt(len(info_scores))) + 10)
    
    counts, bin_edges, patches = ax_main.hist(info_scores, bins=bins, alpha=0.7, 
                                             color='steelblue', edgecolor='black', 
                                             density=True, label=f'Distribution (n={len(info_scores):,})')
    
    # Add kde overlay
    try:
        kde_x = np.linspace(0, 1, 200)
        kde = scipy_stats.gaussian_kde(info_scores)
        kde_y = kde(kde_x)
        ax_main.plot(kde_x, kde_y, 'r-', linewidth=2, label='Kernel Density')
    except:
        pass
    
    # Add vertical lines for common thresholds
    thresholds = [0.3, 0.5, 0.7, 0.8, 0.9]
    threshold_colors = ['red', 'orange', 'yellow', 'lightgreen', 'green']
    
    for thresh, color in zip(thresholds, threshold_colors):
        count_above = (info_scores >= thresh).sum()
        prop_above = count_above / len(info_scores)
        ax_main.axvline(thresh, color=color, linestyle='--', alpha=0.8, 
                       label=f'{thresh}: {prop_above:.1%} above')
    
    ax_main.set_xlabel('INFO Score', fontsize=12)
    ax_main.set_ylabel('Density', fontsize=12)
    ax_main.set_title(f'INFO Score Distribution\\n{output_prefix}', fontsize=14)
    ax_main.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax_main.grid(True, alpha=0.3)

    # Summary statistics box (top right)
    ax_stats = fig.add_subplot(gs[0, 2])
    ax_stats.axis('off')
    
    stats_text = []
    stats_text.append(f"Total variants: {len(info_scores):,}")
    stats_text.append(f"Mean: {np.mean(info_scores):.4f}")
    stats_text.append(f"Median: {np.median(info_scores):.4f}")
    stats_text.append(f"Std Dev: {np.std(info_scores):.4f}")
    stats_text.append(f"Min: {np.min(info_scores):.4f}")
    stats_text.append(f"Max: {np.max(info_scores):.4f}")
    stats_text.append("")
    stats_text.append("Percentiles:")
    for p in [5, 10, 25, 75, 90, 95]:
        val = np.percentile(info_scores, p)
        stats_text.append(f"  {p}th: {val:.4f}")
    
    ax_stats.text(0.1, 0.9, '\\n'.join(stats_text), transform=ax_stats.transAxes, 
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    # Cumulative distribution (middle left)
    ax_cdf = fig.add_subplot(gs[1, 0])
    
    sorted_scores = np.sort(info_scores)
    p_values = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    
    ax_cdf.plot(sorted_scores, p_values, 'b-', linewidth=2, label='Empirical CDF')
    
    # Add threshold lines
    for thresh, color in zip(thresholds, threshold_colors):
        prop_below = (info_scores <= thresh).sum() / len(info_scores)
        ax_cdf.axvline(thresh, color=color, linestyle='--', alpha=0.6)
        ax_cdf.axhline(prop_below, color=color, linestyle=':', alpha=0.6)
    
    ax_cdf.set_xlabel('INFO Score', fontsize=11)
    ax_cdf.set_ylabel('Cumulative Probability', fontsize=11)
    ax_cdf.set_title('Cumulative Distribution', fontsize=12)
    ax_cdf.grid(True, alpha=0.3)
    ax_cdf.legend()

    # Quality categories pie chart (middle center)
    ax_pie = fig.add_subplot(gs[1, 1])
    
    # Define quality categories
    high_quality = (info_scores >= 0.8).sum()
    med_quality = ((info_scores >= 0.5) & (info_scores < 0.8)).sum()
    low_quality = ((info_scores >= 0.3) & (info_scores < 0.5)).sum()
    very_low_quality = (info_scores < 0.3).sum()
    
    categories = ['High\\n(≥0.8)', 'Medium\\n(0.5-0.8)', 'Low\\n(0.3-0.5)', 'Very Low\\n(<0.3)']
    counts = [high_quality, med_quality, low_quality, very_low_quality]
    colors = ['green', 'yellow', 'orange', 'red']
    
    # Only include non-zero categories
    non_zero_idx = [i for i, c in enumerate(counts) if c > 0]
    if non_zero_idx:
        plot_categories = [categories[i] for i in non_zero_idx]
        plot_counts = [counts[i] for i in non_zero_idx]
        plot_colors = [colors[i] for i in non_zero_idx]
        
        wedges, texts, autotexts = ax_pie.pie(plot_counts, labels=plot_categories, 
                                             autopct='%1.1f%%', colors=plot_colors, 
                                             startangle=90)
        
        for autotext in autotexts:
            autotext.set_color('black')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(9)
    
    ax_pie.set_title('Quality Categories', fontsize=12)

    # Box plot by INFO score ranges (middle right)
    ax_box = fig.add_subplot(gs[1, 2])
    
    # Create bins for box plot
    bin_edges = [0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
    bin_labels = ['0-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0']
    
    bin_data = []
    bin_names = []
    for i in range(len(bin_edges) - 1):
        mask = (info_scores >= bin_edges[i]) & (info_scores < bin_edges[i+1])
        if i == len(bin_edges) - 2:  # Include 1.0 in the last bin
            mask = (info_scores >= bin_edges[i]) & (info_scores <= bin_edges[i+1])
        
        if mask.sum() > 0:
            bin_data.append(info_scores[mask])
            bin_names.append(f'{bin_labels[i]}\\n(n={mask.sum()})')
    
    if bin_data:
        bp = ax_box.boxplot(bin_data, labels=bin_names, patch_artist=True)
        
        # Color the boxes
        for patch, color in zip(bp['boxes'], plt.cm.viridis(np.linspace(0.2, 0.9, len(bin_data)))):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
    
    ax_box.set_ylabel('INFO Score', fontsize=11)
    ax_box.set_title('Distribution by Ranges', fontsize=12)
    ax_box.grid(True, alpha=0.3)
    plt.setp(ax_box.get_xticklabels(), rotation=45, ha='right')

    # MAF vs INFO scatter (if MAF available) (bottom left)
    ax_scatter = fig.add_subplot(gs[2, 0])
    
    # Look for MAF column
    maf_cols = ['MAF', 'EmpMAF', 'ALT_Frq', 'AF']
    maf_col = None
    
    for col in maf_cols:
        if col in clean_data.columns:
            maf_col = col
            break
    
    if maf_col:
        # Clean MAF data
        maf_data = pd.to_numeric(clean_data[maf_col], errors='coerce')
        valid_mask = (maf_data.notna() & (maf_data > 0) & (maf_data <= 0.5))
        
        if valid_mask.sum() > 100:
            plot_maf = maf_data[valid_mask]
            plot_info = clean_data.loc[valid_mask, 'info_score']
            
            # Sample if too many points
            if len(plot_maf) > 10000:
                sample_idx = np.random.choice(len(plot_maf), 10000, replace=False)
                plot_maf = plot_maf.iloc[sample_idx]
                plot_info = plot_info.iloc[sample_idx]
            
            scatter = ax_scatter.scatter(plot_maf, plot_info, alpha=0.5, s=10, c=plot_info, 
                                       cmap='viridis', edgecolors='none')
            
            # Add trend line
            try:
                z = np.polyfit(plot_maf, plot_info, 1)
                p = np.poly1d(z)
                x_trend = np.linspace(plot_maf.min(), plot_maf.max(), 100)
                ax_scatter.plot(x_trend, p(x_trend), "r--", alpha=0.8, linewidth=2)
                
                # Calculate correlation
                correlation = np.corrcoef(plot_maf, plot_info)[0, 1]
                ax_scatter.text(0.05, 0.95, f'r = {correlation:.3f}', 
                              transform=ax_scatter.transAxes,
                              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            except:
                pass
            
            ax_scatter.set_xlabel(f'{maf_col}', fontsize=11)
            ax_scatter.set_ylabel('INFO Score', fontsize=11)
            ax_scatter.set_title('INFO vs MAF', fontsize=12)
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax_scatter)
            cbar.set_label('INFO Score', rotation=270, labelpad=15)
        else:
            ax_scatter.text(0.5, 0.5, 'Insufficient MAF data\\nfor scatter plot', 
                           ha='center', va='center', transform=ax_scatter.transAxes)
            ax_scatter.set_title('INFO vs MAF - No Data', fontsize=12)
    else:
        ax_scatter.text(0.5, 0.5, 'MAF data not available', 
                       ha='center', va='center', transform=ax_scatter.transAxes)
        ax_scatter.set_title('INFO vs MAF - No MAF Data', fontsize=12)

    # Violin plot (bottom center)
    ax_violin = fig.add_subplot(gs[2, 1])
    
    if bin_data:
        parts = ax_violin.violinplot(bin_data, positions=range(1, len(bin_data)+1), 
                                   showmeans=True, showmedians=True)
        
        # Customize violin plot
        for pc, color in zip(parts['bodies'], plt.cm.plasma(np.linspace(0.2, 0.9, len(bin_data)))):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        
        ax_violin.set_xticks(range(1, len(bin_names)+1))
        ax_violin.set_xticklabels(bin_names, rotation=45, ha='right')
    
    ax_violin.set_ylabel('INFO Score', fontsize=11)
    ax_violin.set_title('Distribution Shapes', fontsize=12)
    ax_violin.grid(True, alpha=0.3)

    # Summary table (bottom right)
    ax_table = fig.add_subplot(gs[2, 2])
    ax_table.axis('off')
    
    # Create threshold summary
    threshold_data = []
    for thresh in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        count_above = (info_scores >= thresh).sum()
        prop_above = count_above / len(info_scores)
        threshold_data.append([f'≥{thresh}', f'{count_above:,}', f'{prop_above:.1%}'])
    
    table = ax_table.table(cellText=threshold_data,
                          colLabels=['Threshold', 'Count', 'Proportion'],
                          cellLoc='center',
                          loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    # Style the table
    for i in range(len(threshold_data) + 1):
        for j in range(3):
            cell = table[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#4CAF50')
                cell.set_text_props(weight='bold', color='white')
            else:
                cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
    
    ax_table.set_title('Threshold Summary', fontsize=12, pad=20)

    plt.suptitle(f'INFO Score Distribution Analysis: {output_prefix}', fontsize=16, y=0.98)
    plt.savefig(f'{output_prefix}.info_distribution.pdf', dpi=150, bbox_inches='tight')
    plt.close()

    # Calculate and save detailed statistics
    stats_data = []
    stats_data.append(['total_variants', len(info_scores)])
    stats_data.append(['mean', np.mean(info_scores)])
    stats_data.append(['median', np.median(info_scores)])
    stats_data.append(['std', np.std(info_scores)])
    stats_data.append(['min', np.min(info_scores)])
    stats_data.append(['max', np.max(info_scores)])
    stats_data.append(['q25', np.percentile(info_scores, 25)])
    stats_data.append(['q75', np.percentile(info_scores, 75)])
    stats_data.append(['iqr', np.percentile(info_scores, 75) - np.percentile(info_scores, 25)])
    
    # Percentiles
    for p in [1, 5, 10, 90, 95, 99]:
        stats_data.append([f'p{p}', np.percentile(info_scores, p)])
    
    # Threshold counts
    for thresh in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        count_above = (info_scores >= thresh).sum()
        prop_above = count_above / len(info_scores)
        stats_data.append([f'count_above_{thresh}', count_above])
        stats_data.append([f'prop_above_{thresh}', prop_above])
    
    # Quality categories
    stats_data.append(['high_quality_count', high_quality])
    stats_data.append(['med_quality_count', med_quality])
    stats_data.append(['low_quality_count', low_quality])
    stats_data.append(['very_low_quality_count', very_low_quality])
    stats_data.append(['high_quality_prop', high_quality / len(info_scores)])
    stats_data.append(['med_quality_prop', med_quality / len(info_scores)])
    stats_data.append(['low_quality_prop', low_quality / len(info_scores)])
    stats_data.append(['very_low_quality_prop', very_low_quality / len(info_scores)])

    # Save statistics
    stats_df = pd.DataFrame(stats_data, columns=['Statistic', 'Value'])
    stats_df.to_csv(f'{output_prefix}.info_distribution_stats.tsv', sep='\\t', index=False)

    print(f"INFO score distribution analysis completed. Output saved to {output_prefix}.info_distribution.pdf")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub plot
    touch ${prefix}.info_distribution.pdf
    
    # Create stub statistics
    echo -e "Statistic\\tValue" > ${prefix}.info_distribution_stats.tsv
    echo -e "total_variants\\t100000" >> ${prefix}.info_distribution_stats.tsv
    echo -e "mean\\t0.85" >> ${prefix}.info_distribution_stats.tsv
    echo -e "median\\t0.87" >> ${prefix}.info_distribution_stats.tsv
    echo -e "high_quality_prop\\t0.75" >> ${prefix}.info_distribution_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
        seaborn: 0.12.0
        numpy: 1.24.0
        scipy: 1.10.0
    END_VERSIONS
    """
}