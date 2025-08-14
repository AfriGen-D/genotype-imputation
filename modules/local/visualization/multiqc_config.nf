process GENERATE_MULTIQC_CONFIG {
    tag "$meta.id"
    label 'process_low'
    label 'python_plotting'

    input:
    tuple val(meta), path(input_files)
    val(config_params)

    output:
    tuple val(meta), path("*.multiqc_config.yaml"), emit: config
    tuple val(meta), path("multiqc_report.html"), emit: report
    tuple val(meta), path("multiqc_data/"), emit: data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config = config_params ? config_params : [:]
    """
    #!/usr/bin/env python3

    import yaml
    import json
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import subprocess
    import sys
    import os

    # Input parameters
    output_prefix = '${prefix}'
    input_files = [f for f in '${input_files}'.split() if f != 'null' and Path(f).exists()]
    
    print(f"Processing {len(input_files)} input files for MultiQC")

    # Generate MultiQC configuration
    multiqc_config = {
        'title': 'Genotype Imputation Quality Report',
        'subtitle': f'Analysis for: {output_prefix}',
        'intro_text': 'This report summarizes the quality metrics from genotype imputation analysis.',
        
        'report_header_info': [
            {'Sample': output_prefix},
            {'Generated': '${workflow.start}'},
            {'Pipeline': 'genotype-imputation v${workflow.manifest.version}'}
        ],
        
        'show_analysis_paths': False,
        'show_analysis_time': True,
        
        'custom_logo': None,
        'custom_logo_url': None,
        'custom_logo_title': None,
        
        'module_order': [
            {'custom_content': {'name': 'Imputation Overview'}},
            {'custom_content': {'name': 'Quality Metrics'}},
            {'custom_content': {'name': 'Performance Statistics'}},
            {'custom_content': {'name': 'Visualization Summary'}}
        ],
        
        'custom_data': {
            'imputation_metrics': {
                'file_format': 'tsv',
                'section_name': 'Imputation Quality Metrics',
                'description': 'Key performance indicators for genotype imputation quality',
                'plot_type': 'table',
                'pconfig': {
                    'id': 'imputation_metrics_table',
                    'title': 'Imputation Quality Summary',
                    'min': 0,
                    'max': 1,
                    'scale': 'RdYlGn'
                }
            }
        },
        
        'table_columns_visible': {
            'imputation_metrics': {
                'sample': True,
                'mean_r2': True,
                'mean_info': True,
                'total_variants': True,
                'high_quality_prop': True
            }
        },
        
        'table_columns_placement': {
            'imputation_metrics': {
                'mean_r2': 900,
                'mean_info': 800,
                'total_variants': 700,
                'high_quality_prop': 600
            }
        },
        
        'custom_plot_config': {
            'r2_distribution': {
                'id': 'r2_histogram',
                'section_name': 'R² Distribution',
                'description': 'Distribution of imputation accuracy (R²) scores',
                'plot_type': 'linegraph',
                'pconfig': {
                    'title': 'R² Score Distribution',
                    'ylab': 'Frequency',
                    'xlab': 'R² Score'
                }
            },
            
            'info_distribution': {
                'id': 'info_histogram', 
                'section_name': 'INFO Score Distribution',
                'description': 'Distribution of imputation INFO scores',
                'plot_type': 'linegraph',
                'pconfig': {
                    'title': 'INFO Score Distribution',
                    'ylab': 'Frequency', 
                    'xlab': 'INFO Score'
                }
            },
            
            'maf_accuracy': {
                'id': 'maf_accuracy_plot',
                'section_name': 'Accuracy by MAF',
                'description': 'Imputation accuracy across minor allele frequency bins',
                'plot_type': 'bargraph',
                'pconfig': {
                    'title': 'Mean Accuracy by MAF Category',
                    'ylab': 'Mean R² Score',
                    'xlab': 'MAF Category'
                }
            }
        },
        
        'remove_sections': [
            'software_versions'
        ],
        
        'section_comments': {
            'imputation_metrics': 'Quality metrics calculated from imputation INFO files and comparison with reference data.'
        },
        
        'table_cond_formatting_colours': [
            {'blue': 'rgba(52, 152, 219, 0.9)'},
            {'green': 'rgba(46, 125, 50, 0.9)'},
            {'orange': 'rgba(255, 152, 0, 0.9)'},
            {'red': 'rgba(244, 67, 54, 0.9)'}
        ],
        
        'table_cond_formatting_rules': {
            'mean_r2': [
                {'s_eq': 'pass', 'c': 'green'},
                {'s_eq': 'warn', 'c': 'orange'},
                {'s_eq': 'fail', 'c': 'red'}
            ],
            'mean_info': [
                {'s_eq': 'pass', 'c': 'green'},
                {'s_eq': 'warn', 'c': 'orange'},
                {'s_eq': 'fail', 'c': 'red'}
            ]
        }
    }

    # Write MultiQC config
    config_file = f'{output_prefix}.multiqc_config.yaml'
    with open(config_file, 'w') as f:
        yaml.dump(multiqc_config, f, default_flow_style=False, sort_keys=False)
    
    print(f"MultiQC config written to {config_file}")

    # Process input files to create MultiQC-compatible data
    multiqc_data = {}
    summary_metrics = {}
    
    for input_file in input_files:
        file_path = Path(input_file)
        file_name = file_path.name
        
        print(f"Processing file: {file_name}")
        
        try:
            if file_name.endswith('.json'):
                # Process JSON metrics files
                with open(input_file, 'r') as f:
                    data = json.load(f)
                    if isinstance(data, dict):
                        for key, value in data.items():
                            if isinstance(value, (int, float)):
                                summary_metrics[key] = value
                                
            elif file_name.endswith('.tsv'):
                # Process TSV statistics files
                df = pd.read_csv(input_file, sep='\\t')
                
                if 'Statistic' in df.columns and 'Value' in df.columns:
                    # Convert to dict format for MultiQC
                    for _, row in df.iterrows():
                        stat_name = str(row['Statistic'])
                        try:
                            value = float(row['Value'])
                            summary_metrics[stat_name] = value
                        except (ValueError, TypeError):
                            summary_metrics[stat_name] = str(row['Value'])
                            
                elif len(df.columns) >= 2:
                    # Generic table format
                    table_data = df.to_dict('records')
                    multiqc_data[file_path.stem] = table_data
                    
        except Exception as e:
            print(f"Warning: Could not process {file_name}: {e}")
            continue

    # Create MultiQC-compatible summary table
    if summary_metrics:
        summary_table = {
            output_prefix: summary_metrics
        }
        
        # Write summary metrics for MultiQC
        with open('imputation_metrics_mqc.tsv', 'w') as f:
            # Write header
            headers = ['Sample'] + list(summary_metrics.keys())
            f.write('\\t'.join(headers) + '\\n')
            
            # Write data
            values = [output_prefix] + [str(summary_metrics.get(h, '')) for h in headers[1:]]
            f.write('\\t'.join(values) + '\\n')

    # Create custom content files for MultiQC
    
    # 1. Overview section
    overview_content = f'''
    # Imputation Analysis Overview
    
    **Sample ID:** {output_prefix}
    
    **Analysis Summary:**
    '''
    
    if 'total_variants' in summary_metrics:
        overview_content += f"\\n- **Total Variants Analyzed:** {summary_metrics['total_variants']:,}"
    if 'mean_r2' in summary_metrics:
        overview_content += f"\\n- **Mean R² Score:** {summary_metrics['mean_r2']:.4f}"
    if 'mean_info' in summary_metrics:
        overview_content += f"\\n- **Mean INFO Score:** {summary_metrics['mean_info']:.4f}"
    if 'prop_high_quality' in summary_metrics:
        prop_pct = float(summary_metrics['prop_high_quality']) * 100
        overview_content += f"\\n- **High Quality Variants (R² ≥ 0.8):** {prop_pct:.1f}%"
    
    overview_content += '''
    
    **Quality Thresholds:**
    - **High Quality:** R² ≥ 0.8, INFO ≥ 0.8
    - **Medium Quality:** R² ≥ 0.5, INFO ≥ 0.5  
    - **Low Quality:** R² < 0.5, INFO < 0.5
    
    '''
    
    with open('overview_mqc.html', 'w') as f:
        f.write(f'''
        <div class="alert alert-info">
            <h4>Imputation Quality Assessment</h4>
            {overview_content.replace('\\n', '<br>\\n')}
        </div>
        ''')

    # 2. Create quality metrics summary
    if summary_metrics:
        quality_summary = {}
        
        # Categorize metrics
        r2_metrics = {k: v for k, v in summary_metrics.items() if 'r2' in k.lower() and isinstance(v, (int, float))}
        info_metrics = {k: v for k, v in summary_metrics.items() if 'info' in k.lower() and isinstance(v, (int, float))}
        count_metrics = {k: v for k, v in summary_metrics.items() if 'count' in k.lower() or 'variants' in k.lower()}
        
        quality_html = '<div class="row">'
        
        if r2_metrics:
            quality_html += '''
            <div class="col-md-4">
                <h5>R² Statistics</h5>
                <table class="table table-condensed">
            '''
            for metric, value in r2_metrics.items():
                clean_name = metric.replace('_', ' ').title()
                if isinstance(value, float):
                    value_str = f'{value:.4f}'
                else:
                    value_str = f'{value:,}'
                quality_html += f'<tr><td>{clean_name}</td><td>{value_str}</td></tr>'
            quality_html += '</table></div>'
        
        if info_metrics:
            quality_html += '''
            <div class="col-md-4">
                <h5>INFO Statistics</h5>
                <table class="table table-condensed">
            '''
            for metric, value in info_metrics.items():
                clean_name = metric.replace('_', ' ').title()
                if isinstance(value, float):
                    value_str = f'{value:.4f}'
                else:
                    value_str = f'{value:,}'
                quality_html += f'<tr><td>{clean_name}</td><td>{value_str}</td></tr>'
            quality_html += '</table></div>'
        
        if count_metrics:
            quality_html += '''
            <div class="col-md-4">
                <h5>Variant Counts</h5>
                <table class="table table-condensed">
            '''
            for metric, value in count_metrics.items():
                clean_name = metric.replace('_', ' ').title()
                value_str = f'{value:,}' if isinstance(value, (int, float)) else str(value)
                quality_html += f'<tr><td>{clean_name}</td><td>{value_str}</td></tr>'
            quality_html += '</table></div>'
        
        quality_html += '</div>'
        
        with open('quality_metrics_mqc.html', 'w') as f:
            f.write(quality_html)

    # Run MultiQC
    try:
        # Prepare MultiQC command
        cmd = [
            'multiqc',
            '.',
            '--config', config_file,
            '--force',
            '--filename', 'multiqc_report.html',
            '--title', f'Imputation Report: {output_prefix}',
            '--comment', f'Genotype imputation quality assessment for {output_prefix}',
            '--dirs',
            '--dirs-depth', '1'
        ]
        
        print(f"Running MultiQC command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd='.')
        
        if result.returncode == 0:
            print("MultiQC completed successfully")
            print(result.stdout)
        else:
            print(f"MultiQC failed with return code {result.returncode}")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            
            # Create minimal report if MultiQC fails
            create_minimal_report()
    
    except FileNotFoundError:
        print("MultiQC not found. Creating minimal HTML report instead.")
        create_minimal_report()
    except Exception as e:
        print(f"Error running MultiQC: {e}")
        create_minimal_report()

    def create_minimal_report():
        \"\"\"Create a basic HTML report if MultiQC is not available\"\"\"
        html_content = f'''
        <!DOCTYPE html>
        <html>
        <head>
            <title>Imputation Report - {output_prefix}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .metrics {{ margin: 20px 0; }}
                .metric {{ display: inline-block; margin: 10px; padding: 15px; background: #e7f3ff; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Genotype Imputation Report</h1>
                <h2>Sample: {output_prefix}</h2>
            </div>
            
            <div class="metrics">
                <h3>Summary Metrics</h3>
        '''
        
        for metric, value in summary_metrics.items():
            clean_name = metric.replace('_', ' ').title()
            if isinstance(value, float):
                value_str = f'{value:.4f}'
            elif isinstance(value, int):
                value_str = f'{value:,}'
            else:
                value_str = str(value)
            
            html_content += f'''
                <div class="metric">
                    <strong>{clean_name}:</strong> {value_str}
                </div>
            '''
        
        html_content += '''
            </div>
        </body>
        </html>
        '''
        
        with open('multiqc_report.html', 'w') as f:
            f.write(html_content)
        
        # Create minimal data directory
        os.makedirs('multiqc_data', exist_ok=True)
        with open('multiqc_data/multiqc_general_stats.txt', 'w') as f:
            f.write('Sample\\tMetrics\\n')
            f.write(f'{output_prefix}\\t{len(summary_metrics)}\\n')

    # Ensure output files exist
    if not Path('multiqc_report.html').exists():
        create_minimal_report()
    
    if not Path('multiqc_data').exists():
        os.makedirs('multiqc_data', exist_ok=True)

    print("MultiQC process completed.")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        multiqc: \$(multiqc --version 2>/dev/null || echo 'not_available')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub MultiQC config
    cat > ${prefix}.multiqc_config.yaml << 'EOF'
title: 'Genotype Imputation Quality Report'
subtitle: 'Analysis for: ${prefix}'
intro_text: 'This report summarizes the quality metrics from genotype imputation analysis.'
show_analysis_paths: false
show_analysis_time: true
EOF

    # Create stub MultiQC report
    cat > multiqc_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head><title>MultiQC Report - ${prefix}</title></head>
<body>
    <h1>MultiQC Imputation Report</h1>
    <p>Sample: ${prefix}</p>
    <p>Stub report generated for testing purposes.</p>
</body>
</html>
EOF

    # Create stub data directory
    mkdir -p multiqc_data
    echo "Sample	Metrics" > multiqc_data/multiqc_general_stats.txt
    echo "${prefix}	stub_data" >> multiqc_data/multiqc_general_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        multiqc: 1.14
        pandas: 2.0.0
        pyyaml: 6.0
    END_VERSIONS
    """
}