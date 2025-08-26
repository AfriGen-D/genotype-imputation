process PLOT_MAF_CONCORDANCE {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'

    input:
    tuple val(meta), path(info_file)
    tuple val(meta_mask), path(masked_file, stageAs: 'masked_*')

    output:
    tuple val(meta), path("*.concordance_maf.pdf"), emit: plot
    tuple val(meta), path("*.concordance_stats.tsv"), emit: stats
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

    # Input files
    info_file = '${info_file}'
    masked_file = '${masked_file}' if '${masked_file}' != 'null' else None
    output_prefix = '${prefix}'

    # Read data
    try:
        data = read_file(info_file)
        print(f"Loaded {len(data)} variants from info file")
        print(f"Available columns: {list(data.columns)}")
    except Exception as e:
        print(f"Error reading info file: {e}")
        # Create empty plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f'Error reading data\\n{str(e)}', 
                ha='center', va='center', fontsize=12, transform=ax.transAxes)
        ax.set_title('Concordance vs MAF Analysis - Error', fontsize=14)
        plt.savefig(f'{output_prefix}.concordance_maf.pdf', dpi=150, bbox_inches='tight')
        plt.close()
        
        # Create empty stats
        with open(f'{output_prefix}.concordance_stats.tsv', 'w') as f:
            f.write('Statistic\\tValue\\n')
            f.write('Error\\tFailed to read input file\\n')
        exit(1)

    # Check for concordance data availability
    concordance_data = None
    has_concordance = False

    # Try to read masked file if provided
    if masked_file:
        try:
            masked_data = read_file(masked_file)
            if 'concordance' in masked_data.columns and 'MAF' in masked_data.columns:
                concordance_data = masked_data
                has_concordance = True
                print(f"Using masked concordance data with {len(concordance_data)} variants")
        except Exception as e:
            print(f"Could not read masked file: {e}")

    # Try using empirical R-squared as concordance proxy
    if not has_concordance:
        # Look for various R-squared columns
        rsq_cols = ['EmpRsq', 'Rsq', 'LooRsq', 'R2']
        maf_cols = ['MAF', 'EmpMAF', 'ALT_Frq', 'AF']
        
        rsq_col = None
        maf_col = None
        
        for col in rsq_cols:
            if col in data.columns:
                rsq_col = col
                break
        
        for col in maf_cols:
            if col in data.columns:
                maf_col = col
                break
        
        if rsq_col and maf_col:
            # Filter out missing values and convert to numeric
            clean_data = data[(data[rsq_col] != '-') & (data[rsq_col].notna()) & 
                             (data[maf_col] != '-') & (data[maf_col].notna())].copy()
            
            try:
                clean_data['concordance'] = pd.to_numeric(clean_data[rsq_col], errors='coerce')
                clean_data['MAF'] = pd.to_numeric(clean_data[maf_col], errors='coerce')
                clean_data = clean_data.dropna(subset=['concordance', 'MAF'])
                
                if len(clean_data) > 0:
                    concordance_data = clean_data
                    has_concordance = True
                    print(f"Using {rsq_col} as concordance proxy with {len(concordance_data)} variants")
            except Exception as e:
                print(f"Error processing {rsq_col} and {maf_col}: {e}")

    # Create the plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Imputation Concordance vs MAF Analysis', fontsize=16)

    if has_concordance and len(concordance_data) > 0:
        # Ensure MAF is in [0, 0.5] range
        concordance_data['MAF_adj'] = concordance_data['MAF'].apply(lambda x: min(x, 1-x) if x > 0.5 else x)
        
        # Filter extreme values
        concordance_data = concordance_data[
            (concordance_data['MAF_adj'] > 0) & 
            (concordance_data['MAF_adj'] <= 0.5) &
            (concordance_data['concordance'] >= 0) & 
            (concordance_data['concordance'] <= 1)
        ]

        # Plot 1: Scatter plot with trend
        ax = axes[0, 0]
        
        # Sample data if too many points
        plot_data = concordance_data.copy()
        if len(plot_data) > 10000:
            plot_data = plot_data.sample(n=10000, random_state=42)
        
        ax.scatter(plot_data['MAF_adj'], plot_data['concordance'], 
                  alpha=0.4, s=8, c='steelblue', edgecolors='none')
        
        # Add binned trend line
        maf_bins = np.linspace(0, 0.5, 21)
        bin_centers = []
        bin_means = []
        bin_stds = []
        
        for i in range(len(maf_bins)-1):
            mask = ((concordance_data['MAF_adj'] >= maf_bins[i]) & 
                    (concordance_data['MAF_adj'] < maf_bins[i+1]))
            if mask.sum() > 10:  # At least 10 variants per bin
                bin_centers.append((maf_bins[i] + maf_bins[i+1])/2)
                bin_means.append(concordance_data.loc[mask, 'concordance'].mean())
                bin_stds.append(concordance_data.loc[mask, 'concordance'].std())
        
        if bin_centers:
            ax.plot(bin_centers, bin_means, 'r-', linewidth=2.5, label='Binned mean')
            ax.fill_between(bin_centers, 
                           np.array(bin_means) - np.array(bin_stds),
                           np.array(bin_means) + np.array(bin_stds),
                           alpha=0.2, color='red', label='±1 std')
        
        ax.set_xlabel('Minor Allele Frequency', fontsize=11)
        ax.set_ylabel('Concordance Rate', fontsize=11)
        ax.set_title('Concordance vs MAF (Continuous)', fontsize=12)
        ax.set_xlim(-0.01, 0.51)
        ax.set_ylim(0, 1.05)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 2: Binned bar plot
        ax = axes[0, 1]
        
        # Define MAF bins
        maf_bin_edges = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
        maf_labels = ['Rare\\n(0-0.1%)', 'Very Low\\n(0.1-1%)', 'Low\\n(1-5%)', 
                      'Medium\\n(5-10%)', 'Common\\n(10-20%)', 'Very Common\\n(20-50%)']
        
        concordance_data['MAF_bin'] = pd.cut(concordance_data['MAF_adj'], 
                                           bins=maf_bin_edges, labels=maf_labels)
        
        grouped = concordance_data.groupby('MAF_bin')['concordance'].agg(['mean', 'std', 'count'])
        grouped = grouped.dropna()
        
        if len(grouped) > 0:
            x = range(len(grouped))
            colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(grouped)))
            
            bars = ax.bar(x, grouped['mean'], yerr=grouped['std'], capsize=5, 
                         color=colors, edgecolor='black', alpha=0.7)
            
            # Add value labels
            for i, (bar, (idx, row)) in enumerate(zip(bars, grouped.iterrows())):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                       f'{row["mean"]:.3f}\\n(n={int(row["count"])})',
                       ha='center', va='bottom', fontsize=9)
            
            ax.set_xticks(x)
            ax.set_xticklabels(grouped.index, rotation=45, ha='right')
            ax.set_ylabel('Mean Concordance', fontsize=11)
            ax.set_title('Concordance by MAF Bins', fontsize=12)
            ax.set_ylim(0, 1.05)
            
            # Add reference lines
            ax.axhline(y=0.95, color='g', linestyle='--', alpha=0.7, label='95%')
            ax.axhline(y=0.90, color='orange', linestyle='--', alpha=0.7, label='90%')
            ax.axhline(y=0.80, color='r', linestyle='--', alpha=0.7, label='80%')
            ax.legend(loc='lower right', fontsize=9)

        # Plot 3: Distribution histograms
        ax = axes[1, 0]
        
        # Concordance distribution
        ax.hist(concordance_data['concordance'], bins=50, alpha=0.7, color='skyblue', 
               edgecolor='black', label='Concordance')
        ax.set_xlabel('Concordance Rate', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('Distribution of Concordance Rates', fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 4: MAF distribution
        ax = axes[1, 1]
        
        ax.hist(concordance_data['MAF_adj'], bins=50, alpha=0.7, color='lightcoral', 
               edgecolor='black', label='MAF')
        ax.set_xlabel('Minor Allele Frequency', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('Distribution of Minor Allele Frequencies', fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Calculate and save statistics
        stats_data = []
        stats_data.append(['total_variants', len(concordance_data)])
        stats_data.append(['mean_concordance', concordance_data['concordance'].mean()])
        stats_data.append(['median_concordance', concordance_data['concordance'].median()])
        stats_data.append(['std_concordance', concordance_data['concordance'].std()])
        stats_data.append(['min_concordance', concordance_data['concordance'].min()])
        stats_data.append(['max_concordance', concordance_data['concordance'].max()])
        stats_data.append(['mean_maf', concordance_data['MAF_adj'].mean()])
        stats_data.append(['median_maf', concordance_data['MAF_adj'].median()])
        
        # Quality metrics
        high_qual = (concordance_data['concordance'] >= 0.9).sum()
        med_qual = ((concordance_data['concordance'] >= 0.8) & 
                   (concordance_data['concordance'] < 0.9)).sum()
        low_qual = (concordance_data['concordance'] < 0.8).sum()
        
        stats_data.append(['variants_high_concordance_90', high_qual])
        stats_data.append(['variants_med_concordance_80_90', med_qual])
        stats_data.append(['variants_low_concordance_80', low_qual])
        stats_data.append(['prop_high_concordance_90', high_qual / len(concordance_data)])
        
        # By MAF bin statistics
        if len(grouped) > 0:
            for bin_name, row in grouped.iterrows():
                clean_bin_name = str(bin_name).replace('\\n', '_').replace('(', '').replace(')', '').replace('%', 'pct')
                stats_data.append([f'mean_concordance_{clean_bin_name}', row['mean']])
                stats_data.append([f'count_{clean_bin_name}', int(row['count'])])

    else:
        # No concordance data available
        for ax in axes.flat:
            ax.text(0.5, 0.5, 'Concordance data not available\\n\\n' +
                   'To generate concordance metrics:\\n' +
                   '1. Mask known variants before imputation\\n' +
                   '2. Compare imputed vs true genotypes\\n' +
                   '3. Calculate concordance rates',
                   ha='center', va='center', fontsize=11,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            ax.set_title('No Concordance Data', fontsize=12)
        
        stats_data = [['status', 'no_concordance_data']]

    # Save statistics
    stats_df = pd.DataFrame(stats_data, columns=['Statistic', 'Value'])
    stats_df.to_csv(f'{output_prefix}.concordance_stats.tsv', sep='\\t', index=False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.concordance_maf.pdf', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Concordance vs MAF analysis completed. Output saved to {output_prefix}.concordance_maf.pdf")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub plot
    touch ${prefix}.concordance_maf.pdf
    
    # Create stub statistics
    echo -e "Statistic\\tValue" > ${prefix}.concordance_stats.tsv
    echo -e "total_variants\\t50000" >> ${prefix}.concordance_stats.tsv
    echo -e "mean_concordance\\t0.95" >> ${prefix}.concordance_stats.tsv
    echo -e "median_concordance\\t0.96" >> ${prefix}.concordance_stats.tsv
    echo -e "prop_high_concordance_90\\t0.85" >> ${prefix}.concordance_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
        seaborn: 0.12.0
        numpy: 1.24.0
    END_VERSIONS
    """
}