#!/usr/bin/env python3
"""
Generate an HTML index page for all imputation reports
"""

import os
import json
import argparse
from pathlib import Path
from datetime import datetime
import glob

def scan_report_structure(reports_dir):
    """Scan and categorize the report directory structure"""
    reports_dir = Path(reports_dir)
    
    structure = {
        'pdfs': [],
        'jsons': [],
        'plots': {'chunk': [], 'chromosome': [], 'genome': []},
        'warnings': [],
        'summaries': []
    }
    
    # Walk through directory structure
    for item in reports_dir.rglob('*'):
        if item.is_file():
            rel_path = item.relative_to(reports_dir)
            
            # Categorize by file type and location
            if item.suffix == '.pdf':
                structure['pdfs'].append(rel_path)
                
                # Categorize plots by level
                if 'genome' in str(rel_path):
                    structure['plots']['genome'].append(rel_path)
                elif 'chromosome' in str(rel_path) or 'chr' in item.name:
                    structure['plots']['chromosome'].append(rel_path)
                elif 'chunk' in str(rel_path):
                    structure['plots']['chunk'].append(rel_path)
                    
            elif item.suffix == '.json':
                structure['jsons'].append(rel_path)
                if 'summary' in item.name:
                    structure['summaries'].append(rel_path)
                    
            elif 'warning' in str(rel_path):
                structure['warnings'].append(rel_path)
    
    return structure

def generate_html_index(reports_dir, output_file, dataset_name="ChiPImputation"):
    """Generate an HTML index page for all reports"""
    
    # Scan directory structure
    structure = scan_report_structure(reports_dir)
    
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{dataset_name} - Imputation Reports</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        .header {{
            background: white;
            border-radius: 10px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2d3748;
            margin-bottom: 10px;
            font-size: 2.5em;
        }}
        .subtitle {{
            color: #718096;
            font-size: 1.1em;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 30px;
        }}
        .stat-card {{
            background: #f7fafc;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #2d3748;
        }}
        .stat-label {{
            color: #718096;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        .reports-section {{
            background: white;
            border-radius: 10px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }}
        h2 {{
            color: #2d3748;
            margin-bottom: 20px;
            border-bottom: 2px solid #e2e8f0;
            padding-bottom: 10px;
        }}
        .report-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }}
        .report-card {{
            background: #f7fafc;
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 15px;
            transition: all 0.3s ease;
        }}
        .report-card:hover {{
            transform: translateY(-3px);
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            background: #fff;
        }}
        .report-title {{
            font-weight: 600;
            color: #2d3748;
            margin-bottom: 8px;
        }}
        .report-meta {{
            font-size: 0.85em;
            color: #718096;
        }}
        .report-link {{
            display: inline-block;
            margin-top: 10px;
            padding: 6px 12px;
            background: #667eea;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            font-size: 0.9em;
            transition: background 0.3s ease;
        }}
        .report-link:hover {{
            background: #5a67d8;
        }}
        .quality-badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            font-weight: 600;
            margin-left: 8px;
        }}
        .quality-high {{ background: #48bb78; color: white; }}
        .quality-medium {{ background: #ed8936; color: white; }}
        .quality-low {{ background: #f56565; color: white; }}
        .footer {{
            text-align: center;
            color: white;
            padding: 20px;
            margin-top: 40px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 {dataset_name} Imputation Reports</h1>
            <p class="subtitle">Generated on {timestamp}</p>
            
            <div class="stats-grid">
                {stats_cards}
            </div>
        </div>
        
        {report_sections}
        
        <div class="footer">
            <p>Generated by ChiPImputation Pipeline | © 2025 H3ABioNet</p>
        </div>
    </div>
</body>
</html>
    """
    
    # Use structured data from scan
    pdf_files = structure['pdfs']
    json_files = structure['jsons']
    
    # Use pre-categorized reports
    chunk_reports = [(Path(p).stem, p) for p in structure['plots']['chunk']]
    chr_reports = [(Path(p).stem, p) for p in structure['plots']['chromosome']]
    genome_reports = [(Path(p).stem, p) for p in structure['plots']['genome']]
    warnings = structure['warnings']
    summaries = structure['summaries']
    
    # Generate statistics cards
    stats_cards = f"""
        <div class="stat-card">
            <div class="stat-value">{len(pdf_files)}</div>
            <div class="stat-label">Total Reports</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(genome_reports)}</div>
            <div class="stat-label">Genome Reports</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(chr_reports)}</div>
            <div class="stat-label">Chromosome Reports</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(chunk_reports)}</div>
            <div class="stat-label">Chunk Reports</div>
        </div>
    """
    
    # Generate report sections
    report_sections = ""
    
    # Genome-wide reports
    if genome_reports:
        report_sections += """
        <div class="reports-section">
            <h2>📊 Genome-Wide Reports</h2>
            <div class="report-grid">
        """
        for name, path in sorted(genome_reports):
            report_sections += f"""
                <div class="report-card">
                    <div class="report-title">🌍 {name}</div>
                    <div class="report-meta">Full genome analysis</div>
                    <a href="{path}" class="report-link">View Report</a>
                </div>
            """
        report_sections += "</div></div>"
    
    # Chromosome reports
    if chr_reports:
        report_sections += """
        <div class="reports-section">
            <h2>🧬 Chromosome-Level Reports</h2>
            <div class="report-grid">
        """
        for name, path in sorted(chr_reports)[:20]:  # Limit display
            chr_num = name.split('_')[1] if '_' in name else 'Unknown'
            report_sections += f"""
                <div class="report-card">
                    <div class="report-title">Chr {chr_num}</div>
                    <div class="report-meta">{name}</div>
                    <a href="{path}" class="report-link">View Report</a>
                </div>
            """
        if len(chr_reports) > 20:
            report_sections += f"""
                <div class="report-card">
                    <div class="report-title">+ {len(chr_reports) - 20} more</div>
                    <div class="report-meta">Additional chromosome reports available</div>
                </div>
            """
        report_sections += "</div></div>"
    
    # Chunk reports summary
    if chunk_reports:
        report_sections += f"""
        <div class="reports-section">
            <h2>📦 Chunk-Level Reports</h2>
            <p style="margin: 20px 0; color: #718096;">
                {len(chunk_reports)} chunk reports generated. These provide detailed metrics for individual genomic regions.
            </p>
            <div class="report-grid">
        """
        # Show first few chunks as examples
        for name, path in sorted(chunk_reports)[:6]:
            report_sections += f"""
                <div class="report-card">
                    <div class="report-title">{name[:30]}...</div>
                    <div class="report-meta">Chunk analysis</div>
                    <a href="{path}" class="report-link">View Report</a>
                </div>
            """
        report_sections += "</div></div>"
    
    # Warnings section (if any)
    if warnings:
        report_sections += """
        <div class="reports-section" style="background: #fff3cd; border-left: 4px solid #ffc107;">
            <h2>⚠️ Warnings & Issues</h2>
            <p style="margin: 20px 0; color: #856404;">
                {warning_count} warnings detected during processing. Review these for potential issues.
            </p>
            <ul style="list-style: none; padding: 0;">
        """.format(warning_count=len(warnings))
        for warning_path in warnings[:10]:  # Show first 10 warnings
            report_sections += f"""
                <li style="margin: 10px 0;">
                    <a href="{warning_path}" style="color: #856404;">
                        📄 {Path(warning_path).name}
                    </a>
                </li>
            """
        if len(warnings) > 10:
            report_sections += f"""
                <li style="margin: 10px 0; color: #856404;">
                    ... and {len(warnings) - 10} more warnings
                </li>
            """
        report_sections += "</ul></div>"
    
    # Generate final HTML
    html_content = html_template.format(
        dataset_name=dataset_name,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        stats_cards=stats_cards,
        report_sections=report_sections
    )
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"Report index generated: {output_file}")
    return len(pdf_files), len(json_files)

def main():
    parser = argparse.ArgumentParser(description='Generate HTML index for imputation reports')
    parser.add_argument('reports_dir', help='Directory containing reports')
    parser.add_argument('--output', default='report_index.html', help='Output HTML file')
    parser.add_argument('--dataset', default='ChiPImputation', help='Dataset name')
    
    args = parser.parse_args()
    
    pdf_count, json_count = generate_html_index(
        args.reports_dir, 
        args.output,
        args.dataset
    )
    
    print(f"Indexed {pdf_count} PDF reports and {json_count} JSON summaries")

if __name__ == "__main__":
    main()