#!/usr/bin/env python3
"""
Generate index.html files for easy navigation of organized results.
"""

import os
import argparse
from pathlib import Path
from datetime import datetime
import json

def generate_html_header(title):
    """Generate HTML header with styling."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            padding: 30px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2d3748;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #4a5568;
            margin-top: 30px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
        }}
        .stat-label {{
            font-size: 0.9em;
            opacity: 0.9;
            margin-top: 5px;
        }}
        .file-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .file-card {{
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 15px;
            transition: all 0.3s;
        }}
        .file-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }}
        .file-card a {{
            color: #667eea;
            text-decoration: none;
            font-weight: 500;
        }}
        .file-card a:hover {{
            text-decoration: underline;
        }}
        .breadcrumb {{
            background: #f7fafc;
            padding: 10px 15px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .breadcrumb a {{
            color: #667eea;
            text-decoration: none;
        }}
        .breadcrumb a:hover {{
            text-decoration: underline;
        }}
        .report-links {{
            background: #f0f4ff;
            border-left: 4px solid #667eea;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .report-links h3 {{
            margin-top: 0;
            color: #2d3748;
        }}
        .report-links a {{
            display: inline-block;
            margin: 5px 10px 5px 0;
            padding: 8px 15px;
            background: white;
            border: 1px solid #667eea;
            border-radius: 20px;
            color: #667eea;
            text-decoration: none;
            transition: all 0.3s;
        }}
        .report-links a:hover {{
            background: #667eea;
            color: white;
        }}
        .timestamp {{
            color: #718096;
            font-size: 0.9em;
            margin-top: 30px;
            text-align: right;
        }}
    </style>
</head>
<body>
    <div class="container">
"""

def generate_html_footer():
    """Generate HTML footer."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"""
        <div class="timestamp">Generated on {timestamp}</div>
    </div>
</body>
</html>"""

def generate_breadcrumb(path, base_path):
    """Generate breadcrumb navigation."""
    parts = Path(path).relative_to(base_path).parts
    breadcrumb = ['<div class="breadcrumb">']
    breadcrumb.append('<a href="/">Home</a>')
    
    for i, part in enumerate(parts):
        if i < len(parts) - 1:
            link_path = '/' + '/'.join(parts[:i+1])
            breadcrumb.append(f' / <a href="{link_path}/index.html">{part}</a>')
        else:
            breadcrumb.append(f' / <strong>{part}</strong>')
    
    breadcrumb.append('</div>')
    return ''.join(breadcrumb)

def collect_metrics(directory):
    """Collect metrics from JSON files in directory."""
    metrics = {
        'total_variants': 0,
        'well_imputed': 0,
        'mean_r2': [],
        'chunks': 0
    }
    
    for json_file in Path(directory).rglob('*.json'):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                if 'total_variants' in data:
                    metrics['total_variants'] += data['total_variants']
                if 'well_imputed_variants' in data:
                    metrics['well_imputed'] += data['well_imputed_variants']
                if 'mean_r2' in data and data['mean_r2']:
                    metrics['mean_r2'].append(data['mean_r2'])
                if 'chunk_id' in data:
                    metrics['chunks'] += 1
        except:
            continue
    
    if metrics['mean_r2']:
        metrics['mean_r2'] = sum(metrics['mean_r2']) / len(metrics['mean_r2'])
    else:
        metrics['mean_r2'] = 0
    
    return metrics

def generate_index(directory, base_path, level='root'):
    """Generate index.html for a directory."""
    index_path = os.path.join(directory, 'index.html')
    dir_name = os.path.basename(directory)
    
    # Determine the level and generate appropriate title
    if level == 'root':
        title = "Imputation Results"
    elif level == 'dataset':
        title = f"Dataset: {dir_name}"
    elif level == 'refpanel':
        title = f"Reference Panel: {dir_name}"
    else:
        title = dir_name.replace('_', ' ').title()
    
    html = [generate_html_header(title)]
    
    # Add breadcrumb if not root
    if level != 'root':
        html.append(generate_breadcrumb(directory, base_path))
    
    html.append(f"<h1>{title}</h1>")
    
    # Collect and display metrics
    metrics = collect_metrics(directory)
    if any([metrics['total_variants'], metrics['well_imputed'], metrics['chunks']]):
        html.append('<div class="stats-grid">')
        
        if metrics['total_variants']:
            html.append(f'''
                <div class="stat-card">
                    <div class="stat-value">{metrics['total_variants']:,}</div>
                    <div class="stat-label">Total Variants</div>
                </div>
            ''')
        
        if metrics['well_imputed']:
            percentage = (metrics['well_imputed'] / metrics['total_variants'] * 100) if metrics['total_variants'] else 0
            html.append(f'''
                <div class="stat-card">
                    <div class="stat-value">{metrics['well_imputed']:,}</div>
                    <div class="stat-label">Well Imputed ({percentage:.1f}%)</div>
                </div>
            ''')
        
        if metrics['mean_r2']:
            html.append(f'''
                <div class="stat-card">
                    <div class="stat-value">{metrics['mean_r2']:.3f}</div>
                    <div class="stat-label">Mean R²</div>
                </div>
            ''')
        
        if metrics['chunks']:
            html.append(f'''
                <div class="stat-card">
                    <div class="stat-value">{metrics['chunks']}</div>
                    <div class="stat-label">Chunks Processed</div>
                </div>
            ''')
        
        html.append('</div>')
    
    # Find and list reports
    report_files = []
    for ext in ['.html', '.pdf']:
        report_files.extend(Path(directory).glob(f'**/final_report{ext}'))
        report_files.extend(Path(directory).glob(f'**/*report{ext}'))
    
    if report_files:
        html.append('<div class="report-links">')
        html.append('<h3>📊 Reports</h3>')
        for report in report_files:
            rel_path = report.relative_to(directory)
            file_type = 'PDF' if str(report).endswith('.pdf') else 'HTML'
            html.append(f'<a href="{rel_path}">{file_type} Report</a>')
        html.append('</div>')
    
    # List subdirectories
    subdirs = [d for d in Path(directory).iterdir() if d.is_dir() and d.name != '.nextflow']
    if subdirs:
        html.append('<h2>📁 Browse Results</h2>')
        html.append('<div class="file-grid">')
        
        for subdir in sorted(subdirs):
            # Count files in subdirectory
            file_count = sum(1 for _ in subdir.rglob('*') if _.is_file())
            
            # Determine icon based on directory name
            if 'chunk' in subdir.name:
                icon = '🧩'
            elif 'chromosome' in subdir.name:
                icon = '🧬'
            elif 'genome' in subdir.name:
                icon = '🌍'
            elif 'report' in subdir.name:
                icon = '📊'
            elif 'plot' in subdir.name:
                icon = '📈'
            else:
                icon = '📁'
            
            html.append(f'''
                <div class="file-card">
                    <h3>{icon} <a href="{subdir.name}/index.html">{subdir.name}</a></h3>
                    <p>{file_count} files</p>
                </div>
            ''')
        
        html.append('</div>')
    
    # List files in current directory
    files = [f for f in Path(directory).iterdir() if f.is_file() and f.name != 'index.html']
    if files:
        html.append('<h2>📄 Files</h2>')
        html.append('<div class="file-grid">')
        
        for file in sorted(files):
            size = os.path.getsize(file) / 1024  # KB
            size_str = f"{size:.1f} KB" if size < 1024 else f"{size/1024:.1f} MB"
            
            # Determine icon based on file type
            if file.suffix == '.pdf':
                icon = '📑'
            elif file.suffix == '.html':
                icon = '🌐'
            elif file.suffix in ['.png', '.jpg', '.jpeg']:
                icon = '🖼️'
            elif file.suffix == '.json':
                icon = '📋'
            else:
                icon = '📄'
            
            html.append(f'''
                <div class="file-card">
                    <div>{icon} <a href="{file.name}">{file.name}</a></div>
                    <div style="color: #718096; font-size: 0.9em;">{size_str}</div>
                </div>
            ''')
        
        html.append('</div>')
    
    html.append(generate_html_footer())
    
    # Write index.html
    with open(index_path, 'w') as f:
        f.write('\n'.join(html))
    
    print(f"Generated index: {index_path}")
    
    # Recursively generate for subdirectories
    for subdir in subdirs:
        if level == 'root':
            next_level = 'dataset'
        elif level == 'dataset':
            next_level = 'refpanel'
        else:
            next_level = 'analysis'
        
        generate_index(subdir, base_path, next_level)

def main():
    parser = argparse.ArgumentParser(description='Generate index.html files for organized results')
    parser.add_argument('--results-dir', default='/scratch3/users/mamana/results_organized',
                       help='Root directory of organized results')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results_dir):
        print(f"Results directory not found: {args.results_dir}")
        return
    
    print(f"Generating index files for: {args.results_dir}")
    generate_index(args.results_dir, args.results_dir, 'root')
    print("Index generation complete!")

if __name__ == "__main__":
    main()