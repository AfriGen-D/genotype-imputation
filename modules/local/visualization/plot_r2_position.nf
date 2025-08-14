process PLOT_R2_BY_POSITION {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'

    input:
    tuple val(meta), path(info_file)
    tuple val(meta_freq), path(target_file)

    output:
    tuple val(meta), path("*.r2_position.pdf"), emit: plot
    tuple val(meta), path("*.r2_position_stats.tsv"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maf_thresh = args.contains('--maf-threshold') ? args.split('--maf-threshold')[1].split(' ')[1] : '0'
    """
    #!/usr/bin/env python3

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import gzip
    from pathlib import Path

    # Set style for professional-looking plots
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")

    def read_file(filename):
        \"\"\"Read potentially gzipped files\"\"\"
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as f:
                return pd.read_csv(f, sep='\\t')
        else:
            return pd.read_csv(filename, sep='\\t')

    # Parse arguments
    info_file = '${info_file}'
    target_file = '${target_file}'
    output_prefix = '${prefix}'
    maf_thresh = float('${maf_thresh}') if '${maf_thresh}' else 0

    try:
        # Read in the files
        print(f"Reading info file: {info_file}")
        info = read_file(info_file)
        print(f"Info file shape: {info.shape}")
        print(f"Info columns: {list(info.columns)}")
        
        print(f"Reading target file: {target_file}")
        freq = read_file(target_file)
        print(f"Frequency file shape: {freq.shape}")
        print(f"Frequency columns: {list(freq.columns)}")
        
        # Identify required columns from info file
        required_info_cols = ['SNP', 'Rsq', 'Genotyped', 'MAF', 'ALT_Frq']
        available_info_cols = []
        
        # Map column names flexibly
        col_mappings = {
            'SNP': ['SNP', 'ID', 'VARIANT_ID', 'MarkerName'],
            'Rsq': ['Rsq', 'R2', 'LooRsq', 'EmpRsq'],
            'Genotyped': ['Genotyped', 'TYPE', 'Status'],
            'MAF': ['MAF', 'EmpMAF', 'MinorAlleleFrequency'],
            'ALT_Frq': ['ALT_Frq', 'AF', 'AlleleFrequency', 'ALT_FREQ']
        }
        
        info_cols_found = {}
        for req_col in required_info_cols:
            for potential_col in col_mappings.get(req_col, [req_col]):
                if potential_col in info.columns:
                    info_cols_found[req_col] = potential_col
                    break
        
        print(f"Info column mappings: {info_cols_found}")
        
        # Identify required columns from frequency file
        freq_required_cols = ['SNP', 'CHR', 'POS']
        freq_cols_found = {}
        
        freq_col_mappings = {
            'SNP': ['SNP', 'ID', 'VARIANT_ID', 'MarkerName'],
            'CHR': ['CHR', 'CHROM', '#CHROM', 'Chromosome'],
            'POS': ['POS', 'BP', 'Position', 'POSITION']
        }
        
        for req_col in freq_required_cols:
            for potential_col in freq_col_mappings.get(req_col, [req_col]):
                if potential_col in freq.columns:
                    freq_cols_found[req_col] = potential_col
                    break
        
        print(f"Frequency column mappings: {freq_cols_found}")
        
        # Check if we have minimum required columns
        if 'SNP' not in info_cols_found or 'Rsq' not in info_cols_found:
            raise ValueError(f"Missing required info columns. Found: {info_cols_found}")
        
        if 'SNP' not in freq_cols_found or 'CHR' not in freq_cols_found or 'POS' not in freq_cols_found:
            raise ValueError(f"Missing required frequency columns. Found: {freq_cols_found}")
        
        # Extract relevant columns with proper names
        info_subset = info[[info_cols_found[col] for col in info_cols_found.keys()]].copy()
        info_subset.columns = list(info_cols_found.keys())
        
        freq_subset = freq[[freq_cols_found[col] for col in freq_cols_found.keys()]].copy()
        freq_subset.columns = list(freq_cols_found.keys())
        
        # Handle SNP ID format conversion if needed
        if '_' in freq_subset['SNP'].iloc[0] if len(freq_subset) > 0 else False:
            # Split SNP column if it's in CHR_POS_REF_ALT format
            snp_parts = freq_subset['SNP'].str.split('_', expand=True)
            if snp_parts.shape[1] >= 4:
                freq_subset['SNP'] = (freq_subset['CHR'].astype(str) + ':' + 
                                    freq_subset['POS'].astype(str) + ':' + 
                                    snp_parts[2] + ':' + snp_parts[3])
        
        # Merge tables
        print("Merging data...")
        full = pd.merge(freq_subset, info_subset, on='SNP', how='inner')
        print(f"After merge: {len(full)} variants")
        
        if len(full) == 0:
            print("Warning: No matching SNPs found between frequency file and info file.")
            print(f"  Frequency file has {len(freq_subset)} SNPs")
            print(f"  Info file has {len(info_subset)} SNPs")
            if len(freq_subset) > 0 and len(info_subset) > 0:
                print(f"  Sample freq SNPs: {freq_subset['SNP'].head(3).tolist()}")
                print(f"  Sample info SNPs: {info_subset['SNP'].head(3).tolist()}")
        
        # Clean and process the data
        full['CHR'] = full['CHR'].astype(str)
        full['POS'] = pd.to_numeric(full['POS'], errors='coerce')
        
        # Clean R-squared values
        if 'Rsq' in full.columns:
            full = full[full['Rsq'] != '-']
            full['Rsq'] = pd.to_numeric(full['Rsq'], errors='coerce')
            full = full[full['Rsq'].notna()]
        
        print(f"After cleaning R-squared: {len(full)} variants")
        
    except Exception as e:
        print(f"Error processing input files: {e}")
        # Create error plot
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, f'Error loading or processing data\\n{str(e)}', 
                ha='center', va='center', fontsize=12, transform=ax.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        ax.set_title('R² vs Position Plot - Error', fontsize=14)
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Imputation accuracy (R-squared)')
        plt.tight_layout()
        plt.savefig(f'{output_prefix}.r2_position.pdf', dpi=150, bbox_inches='tight')
        plt.close()
        
        # Create error stats
        with open(f'{output_prefix}.r2_position_stats.tsv', 'w') as f:
            f.write('Statistic\\tValue\\n')
            f.write(f'Error\\t{str(e)}\\n')
        
        exit(1)

    try:
        # Process chromosome names for display
        full['CHR_display'] = 'Chromosome ' + full['CHR'].astype(str)
        
        # Separate genotyped and imputed variants
        if 'Genotyped' in full.columns:
            imputed = full[full['Genotyped'] == 'Imputed'].copy()
            genotyped = full[full['Genotyped'] == 'Genotyped'].copy()
        else:
            # If no genotype status, assume all are imputed
            imputed = full.copy()
            genotyped = pd.DataFrame()
        
        print(f"Imputed variants: {len(imputed)}")
        print(f"Genotyped variants: {len(genotyped)}")
        
        # Subsample imputed SNPs for visualization (max ~50,000 points)
        if len(imputed) > 50000:
            n_sample = len(imputed) // (len(imputed) // 50000)
            imputed = imputed.iloc[::n_sample, :].copy()
            print(f"Subsampled imputed to {len(imputed)} variants")
        
        # Subsample genotyped SNPs for rug plot (max ~1,000 points)
        if len(genotyped) > 1000:
            n_sample = len(genotyped) // (len(genotyped) // 1000)
            genotyped = genotyped.iloc[::n_sample, :].copy()
            print(f"Subsampled genotyped to {len(genotyped)} variants")
        
        # Apply MAF threshold if specified
        if maf_thresh > 0 and 'MAF' in imputed.columns:
            af_thresh = maf_thresh / 100
            imputed['MAF_num'] = pd.to_numeric(imputed['MAF'], errors='coerce')
            pre_filter_count = len(imputed)
            imputed = imputed[(imputed['MAF_num'] > af_thresh) & (imputed['MAF_num'] != 0)]
            print(f"After MAF filter (>{af_thresh}): {len(imputed)} variants (removed {pre_filter_count - len(imputed)})")
        
        # Categorize MAF levels for coloring
        def categorize_maf(maf):
            try:
                maf_val = float(maf)
                if 0 < maf_val <= 0.001:
                    return "extreme rare (0,0.001]"
                elif 0.001 < maf_val <= 0.01:
                    return "moderate rare (0.001,0.01]"
                elif 0.01 < maf_val <= 0.02:
                    return "rare (0.01,0.02]"
                elif 0.02 < maf_val <= 0.05:
                    return "moderate (0.02,0.05]"
                elif 0.05 < maf_val <= 0.2:
                    return "common (0.05,0.2]"
                elif 0.2 < maf_val <= 0.5:
                    return "extreme common (0.2,0.5]"
                else:
                    return "unknown"
            except (ValueError, TypeError):
                return "unknown"
        
        if 'MAF' in imputed.columns:
            imputed['MAF_category'] = imputed['MAF'].apply(categorize_maf)
        else:
            imputed['MAF_category'] = 'unknown'
        
        # Get unique chromosomes
        chromosomes = sorted(imputed['CHR_display'].unique()) if len(imputed) > 0 else []
        
        print(f"Chromosomes to plot: {chromosomes}")
        
        # Check if there's data to plot
        if len(imputed) == 0 or len(chromosomes) == 0:
            print("Warning: No imputed data found after filtering. Creating empty plot.")
            fig, ax = plt.subplots(1, 1, figsize=(12, 6))
            ax.text(0.5, 0.5, 'No imputed data available\\nafter filtering', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=14,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            ax.set_xlabel('Position [bp]', fontsize=12)
            ax.set_ylabel('Imputation accuracy (R-squared)', fontsize=12)
            ax.set_title('No Data Available', fontsize=14)
            plt.tight_layout()
            plt.savefig(f'{output_prefix}.r2_position.pdf', dpi=150, bbox_inches='tight')
            plt.close()
            
            # Create empty stats
            with open(f'{output_prefix}.r2_position_stats.tsv', 'w') as f:
                f.write('Statistic\\tValue\\n')
                f.write('total_variants\\t0\\n')
                f.write('chromosomes\\t0\\n')
            exit(0)
        
        # Create the plot
        n_chr = len(chromosomes)
        fig_width = max(12, min(20, 4 * n_chr))  # Scale width with number of chromosomes
        fig, axes = plt.subplots(1, n_chr, figsize=(fig_width, 6), sharey=True)
        
        if n_chr == 1:
            axes = [axes]
        
        # Define colors for MAF categories
        maf_colors = {
            "extreme rare (0,0.001]": "#1f77b4",
            "moderate rare (0.001,0.01]": "#ff7f0e", 
            "rare (0.01,0.02]": "#2ca02c",
            "moderate (0.02,0.05]": "#d62728",
            "common (0.05,0.2]": "#9467bd",
            "extreme common (0.2,0.5]": "#8c564b",
            "unknown": "#7f7f7f"
        }
        
        # Plot statistics
        plot_stats = []
        
        # Plot for each chromosome
        for idx, (ax, chr_name) in enumerate(zip(axes, chromosomes)):
            chr_data = imputed[imputed['CHR_display'] == chr_name].copy()
            chr_genotyped = genotyped[genotyped['CHR_display'] == chr_name].copy() if len(genotyped) > 0 else pd.DataFrame()
            
            print(f"Plotting {chr_name}: {len(chr_data)} imputed, {len(chr_genotyped)} genotyped")
            
            # Plot points colored by MAF category
            for category, color in maf_colors.items():
                cat_data = chr_data[chr_data['MAF_category'] == category]
                if not cat_data.empty:
                    ax.scatter(cat_data['POS'], cat_data['Rsq'], 
                              c=color, s=2, alpha=0.6, label=category, edgecolors='none')
            
            # Add horizontal line at R-squared = 0.3 (common quality threshold)
            ax.axhline(y=0.3, color='red', linestyle='--', alpha=0.7, linewidth=1.5, 
                      label='R² = 0.3 threshold')
            
            # Add additional quality thresholds
            ax.axhline(y=0.8, color='green', linestyle=':', alpha=0.5, linewidth=1, 
                      label='R² = 0.8 (high quality)')
            
            # Add rug plot for genotyped SNPs at bottom
            if not chr_genotyped.empty:
                ax.scatter(chr_genotyped['POS'], np.full(len(chr_genotyped), -0.02), 
                          marker='|', s=15, c='black', alpha=0.4, label='Genotyped SNPs')
            
            # Calculate chromosome statistics
            if not chr_data.empty:
                chr_stats = {
                    'chromosome': chr_name,
                    'n_variants': len(chr_data),
                    'mean_r2': chr_data['Rsq'].mean(),
                    'median_r2': chr_data['Rsq'].median(),
                    'min_pos': chr_data['POS'].min(),
                    'max_pos': chr_data['POS'].max(),
                    'high_quality_prop': (chr_data['Rsq'] >= 0.8).mean(),
                    'acceptable_quality_prop': (chr_data['Rsq'] >= 0.3).mean()
                }
                plot_stats.append(chr_stats)
            
            # Formatting
            ax.set_xlabel('Position [bp]', fontsize=11)
            ax.set_title(chr_name, fontsize=12, pad=10)
            ax.set_ylim(-0.05, 1.05)
            ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax.grid(True, alpha=0.3)
            
            # Format x-axis to show positions in a readable way
            if not chr_data.empty:
                pos_range = chr_data['POS'].max() - chr_data['POS'].min()
                if pos_range > 10**6:
                    # Use scientific notation for large ranges
                    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
            
            # Set y-axis label only on leftmost plot
            if idx == 0:
                ax.set_ylabel('Imputation accuracy (R²)', fontsize=11)
        
        # Add legend (only include categories that were actually plotted)
        handles, labels = [], []
        for ax in axes:
            for handle, label in zip(*ax.get_legend_handles_labels()):
                if label not in labels:
                    handles.append(handle)
                    labels.append(label)
        
        if handles:
            # Position legend outside the plot area
            fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), 
                      title='Categories', fontsize=9)
        
        # Add overall title
        plt.suptitle(f'R² vs Genomic Position: {output_prefix}', fontsize=14, y=0.98)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}.r2_position.pdf', dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved successfully to {output_prefix}.r2_position.pdf")
        
        # Save comprehensive statistics
        stats_data = []
        
        # Overall statistics
        if len(imputed) > 0:
            stats_data.extend([
                ['total_variants', len(imputed)],
                ['total_chromosomes', len(chromosomes)],
                ['overall_mean_r2', imputed['Rsq'].mean()],
                ['overall_median_r2', imputed['Rsq'].median()],
                ['overall_std_r2', imputed['Rsq'].std()],
                ['overall_min_r2', imputed['Rsq'].min()],
                ['overall_max_r2', imputed['Rsq'].max()],
                ['prop_r2_above_0.8', (imputed['Rsq'] >= 0.8).mean()],
                ['prop_r2_above_0.5', (imputed['Rsq'] >= 0.5).mean()],
                ['prop_r2_above_0.3', (imputed['Rsq'] >= 0.3).mean()],
                ['maf_threshold_applied', maf_thresh]
            ])
            
            # MAF category statistics
            if 'MAF_category' in imputed.columns:
                for category in maf_colors.keys():
                    cat_data = imputed[imputed['MAF_category'] == category]
                    if not cat_data.empty:
                        clean_cat = category.replace(' ', '_').replace('(', '').replace(']', '').replace(',', '_')
                        stats_data.extend([
                            [f'{clean_cat}_count', len(cat_data)],
                            [f'{clean_cat}_mean_r2', cat_data['Rsq'].mean()],
                            [f'{clean_cat}_median_r2', cat_data['Rsq'].median()]
                        ])
            
            # Per-chromosome statistics
            for chr_stat in plot_stats:
                chr_clean = chr_stat['chromosome'].replace(' ', '_')
                for key, value in chr_stat.items():
                    if key != 'chromosome':
                        stats_data.append([f'{chr_clean}_{key}', value])
        
        # Save statistics to file
        stats_df = pd.DataFrame(stats_data, columns=['Statistic', 'Value'])
        stats_df.to_csv(f'{output_prefix}.r2_position_stats.tsv', sep='\\t', index=False)
        
        print(f"Statistics saved to {output_prefix}.r2_position_stats.tsv")
        
    except Exception as e:
        print(f"Error during plotting: {e}")
        import traceback
        traceback.print_exc()
        
        # Create error plot
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, f'Error during plotting\\n{str(e)}', 
                ha='center', va='center', fontsize=12, transform=ax.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        ax.set_title('R² vs Position Plot - Plotting Error', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'{output_prefix}.r2_position.pdf', dpi=150, bbox_inches='tight')
        plt.close()
        
        # Create error stats
        with open(f'{output_prefix}.r2_position_stats.tsv', 'w') as f:
            f.write('Statistic\\tValue\\n')
            f.write(f'plotting_error\\t{str(e)}\\n')

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
    touch ${prefix}.r2_position.pdf
    
    # Create stub statistics
    echo -e "Statistic\\tValue" > ${prefix}.r2_position_stats.tsv
    echo -e "total_variants\\t75000" >> ${prefix}.r2_position_stats.tsv
    echo -e "total_chromosomes\\t22" >> ${prefix}.r2_position_stats.tsv
    echo -e "overall_mean_r2\\t0.82" >> ${prefix}.r2_position_stats.tsv
    echo -e "prop_r2_above_0.8\\t0.68" >> ${prefix}.r2_position_stats.tsv

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