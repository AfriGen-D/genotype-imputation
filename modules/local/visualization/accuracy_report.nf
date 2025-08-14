process GENERATE_ACCURACY_REPORT {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'

    input:
    tuple val(meta), path(accuracy_file)
    val(group_name)

    output:
    tuple val(meta), path("*.accuracy_by_maf.tsv"), emit: report
    tuple val(meta), path("*.accuracy_by_maf.pdf"), emit: plot
    tuple val(meta), path("*.accuracy_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def group = group_name ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path
    import sys

    # Set up plotting style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")

    def categorize_maf(maf):
        \"\"\"Categorize MAF into bins\"\"\"
        if 0 < maf <= 0.001:
            return "extreme_rare"
        elif 0.001 < maf <= 0.01:
            return "moderate_rare"
        elif 0.01 < maf <= 0.02:
            return "rare"
        elif 0.02 < maf <= 0.05:
            return "moderate"
        elif 0.05 < maf <= 0.2:
            return "common"
        elif 0.2 < maf <= 0.5:
            return "extreme_common"
        else:
            return "unknown"

    def accuracy_by_maf(accuracy_file, output_prefix, group_name):
        \"\"\"Generate accuracy by MAF report and plot\"\"\"
        
        # Initialize data structure
        maf_categories = {
            'extreme_rare': [],
            'moderate_rare': [],
            'rare': [],
            'moderate': [],
            'common': [],
            'extreme_common': []
        }
        
        # Read accuracy file
        try:
            data = pd.read_csv(accuracy_file, sep='\\t')
        except Exception as e:
            print(f"Error reading accuracy file: {e}")
            return False
            
        # Check for required columns
        required_cols = ['MAF', 'LooRsq']  # Looking for Leave-one-out R-squared
        alt_cols = [['MAF', 'Rsq'], ['EmpMAF', 'EmpRsq'], ['ALT_Frq', 'Rsq']]
        
        maf_col = None
        rsq_col = None
        
        # Find appropriate columns
        if all(col in data.columns for col in required_cols):
            maf_col, rsq_col = required_cols
        else:
            for alt_pair in alt_cols:
                if all(col in data.columns for col in alt_pair):
                    maf_col, rsq_col = alt_pair
                    break
        
        if not maf_col or not rsq_col:
            print(f"Could not find required columns. Available columns: {list(data.columns)}")
            return False
            
        print(f"Using MAF column: {maf_col}, R-squared column: {rsq_col}")
        
        # Process data
        total_variants = 0
        processed_variants = 0
        
        for _, row in data.iterrows():
            try:
                maf = float(row[maf_col])
                rsq = float(row[rsq_col])
                
                # Flip MAF if > 0.5 (get minor allele frequency)
                if maf > 0.5:
                    maf = 1 - maf
                    
                total_variants += 1
                
                # Only include variants with positive R-squared
                if rsq > 0:
                    category = categorize_maf(maf)
                    if category in maf_categories:
                        maf_categories[category].append(rsq)
                        processed_variants += 1
                        
            except (ValueError, TypeError):
                continue
        
        print(f"Processed {processed_variants} variants out of {total_variants} total")
        
        # Calculate statistics for each category
        category_labels = ['(0,0.001]', '(0.001,0.01]', '(0.01,0.02]', 
                          '(0.02,0.05]', '(0.05,0.2]', '(0.2,0.5]']
        category_keys = ['extreme_rare', 'moderate_rare', 'rare', 
                        'moderate', 'common', 'extreme_common']
        
        # Generate main report
        report_lines = []
        report_lines.append('\\t'.join([group_name] + category_labels))
        
        # Calculate means
        means = []
        counts = []
        for key in category_keys:
            if len(maf_categories[key]) > 0:
                mean_acc = np.mean(maf_categories[key])
                means.append(mean_acc)
                counts.append(len(maf_categories[key]))
            else:
                means.append(0.0)
                counts.append(0)
        
        # Format means for output
        mean_str = '\\t'.join([f'{m:.3f}' for m in means])
        report_lines.append(f'{group_name}\\t{mean_str}')
        
        # Save main report
        with open(f'{output_prefix}.accuracy_by_maf.tsv', 'w') as f:
            f.write('\\n'.join(report_lines))
        
        # Generate summary report with counts
        summary_lines = []
        summary_lines.append('\\t'.join([group_name] + category_labels + ['TOTAL']))
        
        count_str = '\\t'.join([str(c) for c in counts])
        total_count = sum(counts)
        summary_lines.append(f'{group_name}\\t{count_str}\\t{total_count}')
        
        with open(f'{output_prefix}.accuracy_summary.tsv', 'w') as f:
            f.write('\\n'.join(summary_lines))
        
        # Generate plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Mean accuracy by MAF category
        non_zero_indices = [i for i, c in enumerate(counts) if c > 0]
        if non_zero_indices:
            x_pos = range(len(non_zero_indices))
            plot_means = [means[i] for i in non_zero_indices]
            plot_labels = [category_labels[i] for i in non_zero_indices]
            plot_counts = [counts[i] for i in non_zero_indices]
            
            colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(plot_means)))
            bars = ax1.bar(x_pos, plot_means, color=colors, alpha=0.8, edgecolor='black')
            
            # Add value labels on bars
            for bar, mean_val, count in zip(bars, plot_means, plot_counts):
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                        f'{mean_val:.3f}\\n(n={count})',
                        ha='center', va='bottom', fontsize=9)
            
            ax1.set_xticks(x_pos)
            ax1.set_xticklabels(plot_labels, rotation=45, ha='right')
            ax1.set_ylabel('Mean Accuracy (R²)', fontsize=11)
            ax1.set_title(f'Imputation Accuracy by MAF Category\\n{group_name}', fontsize=12)
            ax1.set_ylim(0, 1.05)
            
            # Add horizontal reference lines
            ax1.axhline(y=0.8, color='green', linestyle='--', alpha=0.7, label='High quality (0.8)')
            ax1.axhline(y=0.5, color='orange', linestyle='--', alpha=0.7, label='Medium quality (0.5)')
            ax1.legend(loc='upper right', fontsize=9)
        
        # Plot 2: Distribution of variant counts
        if any(c > 0 for c in counts):
            colors = plt.cm.Set3(np.linspace(0, 1, len(category_labels)))
            wedges, texts, autotexts = ax2.pie(counts, labels=category_labels, autopct='%1.1f%%',
                                              colors=colors, startangle=90)
            
            # Improve text visibility
            for autotext in autotexts:
                autotext.set_color('black')
                autotext.set_fontsize(9)
                autotext.set_weight('bold')
            
            ax2.set_title(f'Distribution of Variants by MAF\\n{group_name}', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}.accuracy_by_maf.pdf', dpi=150, bbox_inches='tight')
        plt.close()
        
        return True

    # Execute the analysis
    success = accuracy_by_maf('${accuracy_file}', '${prefix}', '${group}')
    
    if not success:
        # Create empty output files if analysis failed
        with open('${prefix}.accuracy_by_maf.tsv', 'w') as f:
            f.write('${group}\\t(0,0.001]\\t(0.001,0.01]\\t(0.01,0.02]\\t(0.02,0.05]\\t(0.05,0.2]\\t(0.2,0.5]\\n')
            f.write('${group}\\t0.000\\t0.000\\t0.000\\t0.000\\t0.000\\t0.000\\n')
        
        with open('${prefix}.accuracy_summary.tsv', 'w') as f:
            f.write('${group}\\t(0,0.001]\\t(0.001,0.01]\\t(0.01,0.02]\\t(0.02,0.05]\\t(0.05,0.2]\\t(0.2,0.5]\\tTOTAL\\n')
            f.write('${group}\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\n')
        
        # Create empty plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'No valid data for accuracy analysis',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        ax.set_title('Accuracy by MAF Analysis - No Data', fontsize=12)
        plt.savefig('${prefix}.accuracy_by_maf.pdf', dpi=150, bbox_inches='tight')
        plt.close()

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def group = group_name ?: "${meta.id}"
    """
    # Create stub accuracy report
    echo -e "${group}\\t(0,0.001]\\t(0.001,0.01]\\t(0.01,0.02]\\t(0.02,0.05]\\t(0.05,0.2]\\t(0.2,0.5]" > ${prefix}.accuracy_by_maf.tsv
    echo -e "${group}\\t0.750\\t0.820\\t0.850\\t0.880\\t0.920\\t0.950" >> ${prefix}.accuracy_by_maf.tsv
    
    # Create stub summary
    echo -e "${group}\\t(0,0.001]\\t(0.001,0.01]\\t(0.01,0.02]\\t(0.02,0.05]\\t(0.05,0.2]\\t(0.2,0.5]\\tTOTAL" > ${prefix}.accuracy_summary.tsv
    echo -e "${group}\\t1500\\t3200\\t4100\\t5800\\t12000\\t18000\\t44600" >> ${prefix}.accuracy_summary.tsv
    
    # Create stub plot
    touch ${prefix}.accuracy_by_maf.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        numpy: 1.24.0
        matplotlib: 3.7.0
        seaborn: 0.12.0
    END_VERSIONS
    """
}