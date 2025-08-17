process AVERAGE_R2 {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/rsquared/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(info_files)
    
    output:
    tuple val(meta), val(ref_name), path("*.average_r2.txt"), emit: average
    tuple val(meta), val(ref_name), path("*.r2_summary.csv"), emit: summary
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def meanr2_out = "${prefix}_${ref_name}.average_r2.txt"
    def summary_out = "${prefix}_${ref_name}.r2_summary.csv"
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import numpy as np
    import glob
    import sys
    
    # Combine all info files
    all_data = []
    file_stats = []
    
    for info_file in glob.glob("*.info*"):
        try:
            df = pd.read_csv(info_file, sep='\\t')
            if 'Rsq' in df.columns or 'R2' in df.columns:
                # Get R² column
                if 'Rsq' in df.columns:
                    r2_col = 'Rsq'
                else:
                    r2_col = 'R2'
                
                # Calculate statistics for this file
                stats = {
                    'file': info_file,
                    'n_variants': len(df),
                    'mean_r2': df[r2_col].mean(),
                    'median_r2': df[r2_col].median(),
                    'std_r2': df[r2_col].std(),
                    'min_r2': df[r2_col].min(),
                    'max_r2': df[r2_col].max(),
                    'q25_r2': df[r2_col].quantile(0.25),
                    'q75_r2': df[r2_col].quantile(0.75)
                }
                
                # Count well-imputed variants
                for threshold in [0.3, 0.5, 0.8]:
                    stats[f'n_r2_ge_{threshold}'] = (df[r2_col] >= threshold).sum()
                    stats[f'pct_r2_ge_{threshold}'] = (df[r2_col] >= threshold).mean() * 100
                
                file_stats.append(stats)
                all_data.append(df[r2_col].values)
        except Exception as e:
            print(f"Warning: Could not process {info_file}: {e}", file=sys.stderr)
            continue
    
    if all_data:
        # Combine all R² values
        all_r2 = np.concatenate(all_data)
        
        # Calculate overall statistics
        overall_stats = {
            'Total variants': len(all_r2),
            'Mean R²': np.mean(all_r2),
            'Median R²': np.median(all_r2),
            'Std R²': np.std(all_r2),
            'Min R²': np.min(all_r2),
            'Max R²': np.max(all_r2),
            'Q25 R²': np.quantile(all_r2, 0.25),
            'Q75 R²': np.quantile(all_r2, 0.75)
        }
        
        # Write average R² file
        with open("${meanr2_out}", 'w') as f:
            f.write(f"Average R² for ${meta.id} - ${ref_name}\\n")
            f.write("=" * 50 + "\\n\\n")
            
            for key, value in overall_stats.items():
                if 'Total' in key:
                    f.write(f"{key}: {value:,}\\n")
                else:
                    f.write(f"{key}: {value:.4f}\\n")
            
            f.write("\\nR² Thresholds:\\n")
            f.write("-" * 30 + "\\n")
            for threshold in [0.3, 0.5, 0.8]:
                n_pass = (all_r2 >= threshold).sum()
                pct_pass = (all_r2 >= threshold).mean() * 100
                f.write(f"R² ≥ {threshold}: {n_pass:,} variants ({pct_pass:.1f}%)\\n")
            
            f.write("\\nMAF Bins (if available):\\n")
            f.write("-" * 30 + "\\n")
            
            # Try to get MAF-specific stats if MAF column exists
            combined_df = pd.concat([pd.read_csv(f, sep='\\t') for f in glob.glob("*.info*") 
                                     if 'MAF' in pd.read_csv(f, sep='\\t', nrows=1).columns or 
                                        'AF' in pd.read_csv(f, sep='\\t', nrows=1).columns], 
                                    ignore_index=True)
            
            if not combined_df.empty:
                # Get R² and MAF columns
                if 'Rsq' in combined_df.columns:
                    combined_df['r2'] = combined_df['Rsq']
                elif 'R2' in combined_df.columns:
                    combined_df['r2'] = combined_df['R2']
                
                if 'MAF' in combined_df.columns:
                    combined_df['maf'] = combined_df['MAF']
                elif 'AF' in combined_df.columns:
                    combined_df['maf'] = combined_df['AF'].apply(lambda x: min(x, 1-x) if pd.notna(x) else x)
                
                if 'maf' in combined_df.columns:
                    maf_bins = [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5)]
                    for low, high in maf_bins:
                        mask = (combined_df['maf'] >= low) & (combined_df['maf'] < high)
                        if mask.any():
                            mean_r2 = combined_df.loc[mask, 'r2'].mean()
                            n_vars = mask.sum()
                            f.write(f"MAF [{low:.2f}-{high:.2f}): mean R²={mean_r2:.3f}, n={n_vars:,}\\n")
        
        # Write summary CSV
        if file_stats:
            summary_df = pd.DataFrame(file_stats)
            summary_df.to_csv("${summary_out}", index=False)
            print(f"Processed {len(file_stats)} info files")
            print(f"Overall mean R²: {overall_stats['Mean R²']:.4f}")
        
    else:
        # No data found
        with open("${meanr2_out}", 'w') as f:
            f.write("No imputation info data found\\n")
        
        pd.DataFrame().to_csv("${summary_out}", index=False)
    
    # Version info
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
        f.write(f'    numpy: {np.__version__}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Average R²: 0.85" > ${prefix}_${ref_name}.average_r2.txt
    echo "file,mean_r2,n_variants" > ${prefix}_${ref_name}.r2_summary.csv
    echo "test.info,0.85,1000" >> ${prefix}_${ref_name}.r2_summary.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        numpy: 1.24.0
    END_VERSIONS
    """
}