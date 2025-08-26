process GENERATE_PERFORMANCE_REPORT {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'

    input:
    tuple val(meta), path(metrics_json)
    tuple val(meta_data), path(r2_data), path(info_data)

    output:
    tuple val(meta), path("*.performance_report.tsv"), emit: report
    tuple val(meta), path("*.performance_summary.txt"), emit: summary  
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import json
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import sys

    # Load metrics data
    with open('${metrics_json}', 'r') as f:
        metrics = json.load(f)

    # Read R2 and INFO data files
    try:
        r2_data = pd.read_csv('${r2_data}', sep='\\t')
        info_data = pd.read_csv('${info_data}', sep='\\t')
        
        # Merge data if both have common identifiers
        if 'SNP' in r2_data.columns and 'SNP' in info_data.columns:
            merged_data = pd.merge(r2_data, info_data, on='SNP', how='inner')
        else:
            merged_data = pd.concat([r2_data, info_data], axis=1)
            
    except Exception as e:
        print(f"Warning: Could not read data files: {e}")
        merged_data = pd.DataFrame()

    # Generate performance report
    report_data = []
    
    # Basic metrics from JSON
    if 'r2_mean' in metrics:
        report_data.append(['R2_mean', metrics['r2_mean'], 'Mean R-squared across all variants'])
    if 'r2_median' in metrics:
        report_data.append(['R2_median', metrics['r2_median'], 'Median R-squared across all variants'])
    if 'info_mean' in metrics:
        report_data.append(['INFO_mean', metrics['info_mean'], 'Mean INFO score across all variants'])
    if 'info_median' in metrics:
        report_data.append(['INFO_median', metrics['info_median'], 'Median INFO score across all variants'])

    # Calculate additional metrics from data if available
    if not merged_data.empty:
        if 'R2' in merged_data.columns or 'Rsq' in merged_data.columns:
            r2_col = 'R2' if 'R2' in merged_data.columns else 'Rsq'
            r2_values = pd.to_numeric(merged_data[r2_col], errors='coerce')
            r2_values = r2_values.dropna()
            
            if len(r2_values) > 0:
                report_data.append(['R2_std', r2_values.std(), 'Standard deviation of R-squared'])
                report_data.append(['R2_min', r2_values.min(), 'Minimum R-squared'])
                report_data.append(['R2_max', r2_values.max(), 'Maximum R-squared'])
                report_data.append(['R2_q25', r2_values.quantile(0.25), '25th percentile R-squared'])
                report_data.append(['R2_q75', r2_values.quantile(0.75), '75th percentile R-squared'])
                
                # Quality thresholds
                r2_high_quality = (r2_values >= 0.8).sum()
                r2_medium_quality = ((r2_values >= 0.5) & (r2_values < 0.8)).sum()
                r2_low_quality = (r2_values < 0.5).sum()
                total_variants = len(r2_values)
                
                report_data.append(['variants_total', total_variants, 'Total number of variants'])
                report_data.append(['variants_high_r2', r2_high_quality, 'Variants with R2 >= 0.8'])
                report_data.append(['variants_medium_r2', r2_medium_quality, 'Variants with 0.5 <= R2 < 0.8'])
                report_data.append(['variants_low_r2', r2_low_quality, 'Variants with R2 < 0.5'])
                report_data.append(['prop_high_r2', r2_high_quality/total_variants, 'Proportion with R2 >= 0.8'])
                report_data.append(['prop_medium_r2', r2_medium_quality/total_variants, 'Proportion with 0.5 <= R2 < 0.8'])
                report_data.append(['prop_low_r2', r2_low_quality/total_variants, 'Proportion with R2 < 0.5'])

        if 'INFO' in merged_data.columns or 'AvgCall' in merged_data.columns:
            info_col = 'INFO' if 'INFO' in merged_data.columns else 'AvgCall'
            info_values = pd.to_numeric(merged_data[info_col], errors='coerce')
            info_values = info_values.dropna()
            
            if len(info_values) > 0:
                report_data.append(['INFO_std', info_values.std(), 'Standard deviation of INFO scores'])
                report_data.append(['INFO_min', info_values.min(), 'Minimum INFO score'])
                report_data.append(['INFO_max', info_values.max(), 'Maximum INFO score'])
                report_data.append(['INFO_q25', info_values.quantile(0.25), '25th percentile INFO score'])
                report_data.append(['INFO_q75', info_values.quantile(0.75), '75th percentile INFO score'])
                
                # INFO quality thresholds
                info_high = (info_values >= 0.8).sum()
                info_medium = ((info_values >= 0.4) & (info_values < 0.8)).sum()
                info_low = (info_values < 0.4).sum()
                
                report_data.append(['variants_high_info', info_high, 'Variants with INFO >= 0.8'])
                report_data.append(['variants_medium_info', info_medium, 'Variants with 0.4 <= INFO < 0.8'])
                report_data.append(['variants_low_info', info_low, 'Variants with INFO < 0.4'])

    # Create performance report DataFrame
    performance_df = pd.DataFrame(report_data, columns=['Metric', 'Value', 'Description'])
    
    # Save performance report
    performance_df.to_csv('${prefix}.performance_report.tsv', sep='\\t', index=False)

    # Generate summary text
    summary_lines = []
    summary_lines.append(f"Imputation Performance Summary for {meta['id'] if isinstance(meta, dict) else meta}")
    summary_lines.append("=" * 60)
    summary_lines.append("")
    
    # Extract key metrics for summary
    key_metrics = ['R2_mean', 'R2_median', 'INFO_mean', 'INFO_median', 
                   'variants_total', 'prop_high_r2', 'prop_high_info']
    
    for _, row in performance_df.iterrows():
        if row['Metric'] in key_metrics:
            if isinstance(row['Value'], float):
                value_str = f"{row['Value']:.4f}"
            else:
                value_str = str(row['Value'])
            summary_lines.append(f"{row['Description']}: {value_str}")
    
    summary_lines.append("")
    summary_lines.append("Quality Classification:")
    summary_lines.append("- High quality (R2 >= 0.8): Well-imputed variants")
    summary_lines.append("- Medium quality (0.5 <= R2 < 0.8): Moderately-imputed variants")
    summary_lines.append("- Low quality (R2 < 0.5): Poorly-imputed variants")
    
    with open('${prefix}.performance_summary.txt', 'w') as f:
        f.write('\\n'.join(summary_lines))

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub performance report
    echo -e "Metric\\tValue\\tDescription" > ${prefix}.performance_report.tsv
    echo -e "R2_mean\\t0.85\\tMean R-squared across all variants" >> ${prefix}.performance_report.tsv
    echo -e "INFO_mean\\t0.90\\tMean INFO score across all variants" >> ${prefix}.performance_report.tsv
    echo -e "variants_total\\t100000\\tTotal number of variants" >> ${prefix}.performance_report.tsv
    
    # Create stub summary
    echo "Imputation Performance Summary" > ${prefix}.performance_summary.txt
    echo "Mean R-squared across all variants: 0.85" >> ${prefix}.performance_summary.txt
    echo "Mean INFO score across all variants: 0.90" >> ${prefix}.performance_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        numpy: 1.24.0
    END_VERSIONS
    """
}