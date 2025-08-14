process CREATE_HTML_REPORT {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'

    input:
    tuple val(meta), path(metrics_files)
    tuple val(meta_plots), path(plot_files)
    tuple val(meta_stats), path(stats_files)
    val(report_title)

    output:
    tuple val(meta), path("*.html"), emit: html_report
    tuple val(meta), path("*.json"), emit: report_data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def title = report_title ?: "Imputation Quality Report"
    """
    #!/usr/bin/env python3

    import json
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import base64
    from datetime import datetime
    import sys

    def encode_image_to_base64(image_path):
        \"\"\"Encode image file to base64 for embedding in HTML\"\"\"
        try:
            with open(image_path, 'rb') as img_file:
                return base64.b64encode(img_file.read()).decode('utf-8')
        except Exception as e:
            print(f"Warning: Could not encode image {image_path}: {e}")
            return None

    def read_tsv_file(file_path):
        \"\"\"Read TSV file with error handling\"\"\"
        try:
            return pd.read_csv(file_path, sep='\\t')
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
            return pd.DataFrame()

    def read_json_file(file_path):
        \"\"\"Read JSON file with error handling\"\"\"
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
            return {}

    # Input parameters
    report_title = '${title}'
    output_prefix = '${prefix}'
    
    # Collect all input files
    all_files = []
    
    # Add metrics files
    metrics_files = [f for f in '${metrics_files}'.split() if f != 'null' and Path(f).exists()]
    all_files.extend(metrics_files)
    
    # Add plot files  
    plot_files = [f for f in '${plot_files}'.split() if f != 'null' and Path(f).exists()]
    all_files.extend(plot_files)
    
    # Add stats files
    stats_files = [f for f in '${stats_files}'.split() if f != 'null' and Path(f).exists()]
    all_files.extend(stats_files)
    
    print(f"Processing {len(all_files)} input files")
    print(f"Metrics files: {len(metrics_files)}")
    print(f"Plot files: {len(plot_files)}")
    print(f"Stats files: {len(stats_files)}")

    # Organize files by type
    pdf_files = [f for f in all_files if f.endswith('.pdf')]
    json_files = [f for f in all_files if f.endswith('.json')]
    tsv_files = [f for f in all_files if f.endswith('.tsv')]
    
    print(f"Found {len(pdf_files)} PDF files, {len(json_files)} JSON files, {len(tsv_files)} TSV files")

    # Read all data
    all_data = {
        'plots': {},
        'metrics': {},
        'statistics': {},
        'summary': {}
    }

    # Process PDF plots (encode to base64 for embedding)
    for pdf_file in pdf_files:
        plot_name = Path(pdf_file).stem
        encoded_plot = encode_image_to_base64(pdf_file)
        if encoded_plot:
            all_data['plots'][plot_name] = {
                'filename': pdf_file,
                'base64': encoded_plot,
                'type': 'pdf'
            }

    # Process JSON metrics
    for json_file in json_files:
        metrics_name = Path(json_file).stem
        metrics_data = read_json_file(json_file)
        if metrics_data:
            all_data['metrics'][metrics_name] = metrics_data

    # Process TSV statistics
    for tsv_file in tsv_files:
        stats_name = Path(tsv_file).stem
        stats_data = read_tsv_file(tsv_file)
        if not stats_data.empty:
            # Convert to dict for JSON serialization
            all_data['statistics'][stats_name] = stats_data.to_dict('records')

    # Generate summary statistics
    summary_stats = {}
    
    # Extract key metrics from different sources
    for metrics_name, metrics in all_data['metrics'].items():
        if 'r2_mean' in metrics:
            summary_stats['mean_r2'] = metrics['r2_mean']
        if 'info_mean' in metrics:
            summary_stats['mean_info'] = metrics['info_mean']
    
    # Extract from statistics files
    for stats_name, stats_list in all_data['statistics'].items():
        if 'performance' in stats_name.lower():
            for stat in stats_list:
                if stat.get('Metric') == 'variants_total':
                    summary_stats['total_variants'] = stat.get('Value', 0)
                elif stat.get('Metric') == 'prop_high_r2':
                    summary_stats['prop_high_quality'] = stat.get('Value', 0)
    
    all_data['summary'] = summary_stats

    # Save comprehensive data as JSON
    with open(f'{output_prefix}.json', 'w') as f:
        json.dump(all_data, f, indent=2, default=str)

    # Generate HTML report
    html_template = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{report_title}</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f5f5f5;
                color: #333;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 30px;
                border-radius: 10px;
                box-shadow: 0 0 20px rgba(0,0,0,0.1);
            }}
            .header {{
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 20px;
                border-bottom: 3px solid #4CAF50;
            }}
            .header h1 {{
                color: #2E7D32;
                margin-bottom: 10px;
                font-size: 2.5em;
            }}
            .header .subtitle {{
                color: #666;
                font-size: 1.2em;
            }}
            .summary {{
                background: linear-gradient(135deg, #E3F2FD 0%, #BBDEFB 100%);
                padding: 25px;
                border-radius: 10px;
                margin-bottom: 30px;
                border-left: 5px solid #2196F3;
            }}
            .summary h2 {{
                color: #1565C0;
                margin-top: 0;
                font-size: 1.8em;
            }}
            .summary-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-top: 20px;
            }}
            .summary-item {{
                background: white;
                padding: 15px;
                border-radius: 8px;
                text-align: center;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            }}
            .summary-item .value {{
                font-size: 2em;
                font-weight: bold;
                color: #2E7D32;
                display: block;
            }}
            .summary-item .label {{
                color: #666;
                font-size: 0.9em;
                margin-top: 5px;
            }}
            .section {{
                margin-bottom: 40px;
            }}
            .section h2 {{
                color: #2E7D32;
                border-bottom: 2px solid #4CAF50;
                padding-bottom: 10px;
                font-size: 1.8em;
            }}
            .section h3 {{
                color: #1565C0;
                margin-top: 25px;
                font-size: 1.4em;
            }}
            .plot-container {{
                text-align: center;
                margin: 25px 0;
                padding: 20px;
                background-color: #f9f9f9;
                border-radius: 10px;
                border: 1px solid #e0e0e0;
            }}
            .plot-container img {{
                max-width: 100%;
                height: auto;
                border-radius: 5px;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            }}
            .plot-title {{
                font-weight: bold;
                margin-bottom: 10px;
                color: #2E7D32;
                font-size: 1.2em;
            }}
            .statistics-table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                background-color: white;
                border-radius: 8px;
                overflow: hidden;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            .statistics-table th {{
                background: linear-gradient(135deg, #4CAF50 0%, #45a049 100%);
                color: white;
                padding: 15px;
                text-align: left;
                font-weight: bold;
            }}
            .statistics-table td {{
                padding: 12px 15px;
                border-bottom: 1px solid #eee;
            }}
            .statistics-table tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            .statistics-table tr:hover {{
                background-color: #f0f8f0;
            }}
            .metrics-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }}
            .metric-card {{
                background: white;
                border: 1px solid #e0e0e0;
                border-radius: 10px;
                padding: 20px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.05);
            }}
            .metric-card h4 {{
                color: #2E7D32;
                margin-top: 0;
                font-size: 1.2em;
            }}
            .footer {{
                text-align: center;
                margin-top: 50px;
                padding-top: 20px;
                border-top: 1px solid #e0e0e0;
                color: #666;
                font-size: 0.9em;
            }}
            .timestamp {{
                color: #999;
                font-size: 0.8em;
                margin-top: 10px;
            }}
            .collapsible {{
                background-color: #f1f1f1;
                color: #333;
                cursor: pointer;
                padding: 15px;
                width: 100%;
                border: none;
                text-align: left;
                outline: none;
                font-size: 1.1em;
                border-radius: 5px;
                margin: 5px 0;
                transition: background-color 0.3s;
            }}
            .collapsible:hover {{
                background-color: #ddd;
            }}
            .collapsible.active {{
                background-color: #4CAF50;
                color: white;
            }}
            .content {{
                padding: 0 15px;
                display: none;
                overflow: hidden;
                background-color: #f9f9f9;
                border-radius: 0 0 5px 5px;
            }}
            .content.show {{
                display: block;
                padding: 15px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>{report_title}</h1>
                <div class="subtitle">Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
            </div>
    '''

    # Add summary section
    if summary_stats:
        html_template += '''
            <div class="summary">
                <h2>📊 Executive Summary</h2>
                <div class="summary-grid">
        '''
        
        if 'total_variants' in summary_stats:
            html_template += f'''
                    <div class="summary-item">
                        <span class="value">{summary_stats['total_variants']:,}</span>
                        <div class="label">Total Variants</div>
                    </div>
            '''
        
        if 'mean_r2' in summary_stats:
            html_template += f'''
                    <div class="summary-item">
                        <span class="value">{summary_stats['mean_r2']:.3f}</span>
                        <div class="label">Mean R²</div>
                    </div>
            '''
        
        if 'mean_info' in summary_stats:
            html_template += f'''
                    <div class="summary-item">
                        <span class="value">{summary_stats['mean_info']:.3f}</span>
                        <div class="label">Mean INFO Score</div>
                    </div>
            '''
        
        if 'prop_high_quality' in summary_stats:
            prop_pct = float(summary_stats['prop_high_quality']) * 100
            html_template += f'''
                    <div class="summary-item">
                        <span class="value">{prop_pct:.1f}%</span>
                        <div class="label">High Quality (R² ≥ 0.8)</div>
                    </div>
            '''
        
        html_template += '''
                </div>
            </div>
        '''

    # Add plots section
    if all_data['plots']:
        html_template += '''
            <div class="section">
                <h2>📈 Visualization Results</h2>
        '''
        
        plot_categories = {
            'distribution': '📊 Score Distributions',
            'position': '🗺️ Genomic Position Analysis',
            'concordance': '🎯 Concordance Analysis',
            'performance': '⚡ Performance Metrics',
            'accuracy': '🔍 Accuracy Analysis'
        }
        
        for category, title in plot_categories.items():
            category_plots = {k: v for k, v in all_data['plots'].items() 
                             if category in k.lower()}
            
            if category_plots:
                html_template += f'''
                <button class="collapsible">{title}</button>
                <div class="content">
                '''
                
                for plot_name, plot_data in category_plots.items():
                    clean_name = plot_name.replace('_', ' ').title()
                    html_template += f'''
                    <div class="plot-container">
                        <div class="plot-title">{clean_name}</div>
                        <embed src="data:application/pdf;base64,{plot_data['base64']}" 
                               type="application/pdf" width="100%" height="600px" />
                        <p style="font-size: 0.9em; color: #666; margin-top: 10px;">
                            <em>File: {plot_data['filename']}</em>
                        </p>
                    </div>
                    '''
                
                html_template += '</div>'
        
        html_template += '</div>'

    # Add statistics section
    if all_data['statistics']:
        html_template += '''
            <div class="section">
                <h2>📋 Detailed Statistics</h2>
        '''
        
        for stats_name, stats_data in all_data['statistics'].items():
            if stats_data:
                clean_name = stats_name.replace('_', ' ').title()
                html_template += f'''
                <button class="collapsible">📊 {clean_name}</button>
                <div class="content">
                    <table class="statistics-table">
                        <thead>
                            <tr>
                '''
                
                # Get column headers
                if isinstance(stats_data[0], dict):
                    headers = list(stats_data[0].keys())
                    for header in headers:
                        html_template += f'<th>{header}</th>'
                    
                    html_template += '''
                            </tr>
                        </thead>
                        <tbody>
                    '''
                    
                    # Add data rows
                    for row in stats_data:
                        html_template += '<tr>'
                        for header in headers:
                            value = row.get(header, '')
                            # Format numeric values
                            if isinstance(value, (int, float)):
                                if isinstance(value, float) and abs(value) < 1:
                                    value = f'{value:.4f}'
                                elif isinstance(value, float):
                                    value = f'{value:.2f}'
                                else:
                                    value = f'{value:,}'
                            html_template += f'<td>{value}</td>'
                        html_template += '</tr>'
                    
                    html_template += '''
                        </tbody>
                    </table>
                </div>
                '''

        html_template += '</div>'

    # Add metrics section
    if all_data['metrics']:
        html_template += '''
            <div class="section">
                <h2>🔬 Metrics Details</h2>
                <div class="metrics-grid">
        '''
        
        for metrics_name, metrics_data in all_data['metrics'].items():
            clean_name = metrics_name.replace('_', ' ').title()
            html_template += f'''
                <div class="metric-card">
                    <h4>{clean_name}</h4>
                    <ul style="list-style-type: none; padding: 0;">
            '''
            
            for key, value in metrics_data.items():
                if isinstance(value, (int, float)):
                    if isinstance(value, float):
                        value = f'{value:.4f}'
                    else:
                        value = f'{value:,}'
                formatted_key = key.replace('_', ' ').title()
                html_template += f'<li><strong>{formatted_key}:</strong> {value}</li>'
            
            html_template += '''
                    </ul>
                </div>
            '''
        
        html_template += '</div></div>'

    # Add footer
    html_template += f'''
            <div class="footer">
                <p>Report generated by Nextflow Genotype Imputation Pipeline</p>
                <div class="timestamp">Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}</div>
            </div>
        </div>

        <script>
            // Add interactivity for collapsible sections
            var coll = document.getElementsByClassName("collapsible");
            var i;

            for (i = 0; i < coll.length; i++) {{
                coll[i].addEventListener("click", function() {{
                    this.classList.toggle("active");
                    var content = this.nextElementSibling;
                    if (content.classList.contains("show")) {{
                        content.classList.remove("show");
                    }} else {{
                        content.classList.add("show");
                    }}
                }});
            }}
        </script>
    </body>
    </html>
    '''

    # Write HTML report
    with open(f'{output_prefix}.html', 'w') as f:
        f.write(html_template)

    print(f"HTML report generated: {output_prefix}.html")
    print(f"Report data saved: {output_prefix}.json")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def title = report_title ?: "Imputation Quality Report"
    """
    # Create stub HTML report
    cat > ${prefix}.html << 'EOF'
<!DOCTYPE html>
<html>
<head><title>${title}</title></head>
<body>
    <h1>${title}</h1>
    <p>Stub HTML report generated for testing purposes.</p>
    <p>Sample ID: ${meta.id}</p>
    <p>Generated: \$(date)</p>
</body>
</html>
EOF

    # Create stub JSON data
    cat > ${prefix}.json << 'EOF'
{
    "plots": {},
    "metrics": {"r2_mean": 0.85, "info_mean": 0.90},
    "statistics": {},
    "summary": {"total_variants": 100000, "mean_r2": 0.85}
}
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
    END_VERSIONS
    """
}