#!/usr/bin/env python3
"""
Generate comprehensive HTML and PDF reports for imputation results.
Combines genome-wide summary data and plots into formatted reports.
"""

import json
import argparse
from pathlib import Path
from datetime import datetime
import base64
import logging
from io import BytesIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings

# Try to import pdfkit, but make it optional
try:
    import pdfkit
    PDFKIT_AVAILABLE = True
except ImportError:
    PDFKIT_AVAILABLE = False
    pdfkit = None

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning)


def load_genome_summary(summary_path):
    """Load genome summary JSON."""
    with open(summary_path, 'r') as f:
        return json.load(f)


def encode_pdf_as_base64(pdf_path):
    """Encode PDF file as base64 for embedding in HTML."""
    try:
        with open(pdf_path, 'rb') as f:
            return base64.b64encode(f.read()).decode('utf-8')
    except Exception as e:
        logger.warning(f"Could not encode PDF {pdf_path}: {e}")
        return None


def format_number(value):
    """Format numbers for display."""
    if isinstance(value, int):
        return f"{value:,}"
    elif isinstance(value, float):
        if value < 0.01:
            return f"{value:.2e}"
        elif value < 1:
            return f"{value:.4f}"
        else:
            return f"{value:,.2f}"
    return str(value)


def generate_html_report(genome_summary, genome_plots_path, dataset, ref_name, population, study):
    """Generate HTML report content."""
    
    # Get current timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Start HTML document
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Imputation Report - {dataset}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
        }}
        .header .subtitle {{
            opacity: 0.9;
            margin-top: 10px;
        }}
        .section {{
            background: white;
            padding: 25px;
            margin-bottom: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .metric-card {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .metric-card .label {{
            font-size: 0.9em;
            color: #666;
            margin-bottom: 5px;
        }}
        .metric-card .value {{
            font-size: 1.5em;
            font-weight: bold;
            color: #333;
        }}
        .table-responsive {{
            overflow-x: auto;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #667eea;
            color: white;
            font-weight: 600;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .good {{
            color: #28a745;
            font-weight: bold;
        }}
        .warning {{
            color: #ffc107;
            font-weight: bold;
        }}
        .poor {{
            color: #dc3545;
            font-weight: bold;
        }}
        .plot-container {{
            text-align: center;
            margin: 30px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            border-top: 1px solid #ddd;
            margin-top: 50px;
        }}
        .chromosome-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .chromosome-card {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            border: 1px solid #dee2e6;
        }}
        .chromosome-card h4 {{
            margin: 0 0 10px 0;
            color: #495057;
        }}
        .progress-bar {{
            width: 100%;
            height: 20px;
            background-color: #e9ecef;
            border-radius: 10px;
            overflow: hidden;
            margin: 10px 0;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea, #764ba2);
            transition: width 0.3s ease;
        }}
        .summary-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        @media print {{
            .header {{
                background: #667eea !important;
                -webkit-print-color-adjust: exact;
                print-color-adjust: exact;
            }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Genotype Imputation Report</h1>
        <div class="subtitle">
            <strong>Dataset:</strong> {dataset} | 
            <strong>Reference:</strong> {ref_name} | 
            <strong>Generated:</strong> {timestamp}
        </div>
    </div>
"""
    
    # Executive Summary Section
    html += f"""
    <div class="section">
        <h2>Executive Summary</h2>
        <div class="summary-box">
            <p>This report summarizes the results of genotype imputation performed on dataset <strong>{dataset}</strong> 
            using reference panel <strong>{ref_name}</strong>.</p>
            <p><strong>Study:</strong> {study}</p>
            <p><strong>Population:</strong> {population}</p>
        </div>
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="label">Chromosomes Processed</div>
                <div class="value">{genome_summary.get('chromosomes_processed', 0)}</div>
            </div>
            <div class="metric-card">
                <div class="label">Total Variants</div>
                <div class="value">{format_number(genome_summary.get('total_variants', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Well-Imputed Variants</div>
                <div class="value">{format_number(genome_summary.get('well_imputed_variants', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Mean R²</div>
                <div class="value {get_quality_class(genome_summary.get('mean_r2', 0))}">{format_number(genome_summary.get('mean_r2', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Mean Concordance</div>
                <div class="value {get_quality_class(genome_summary.get('mean_concordance', 0))}">{format_number(genome_summary.get('mean_concordance', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Imputation Rate</div>
                <div class="value">{format_number((genome_summary.get('imputation_rate', 0) or 0) * 100)}%</div>
            </div>
        </div>
    </div>
"""
    
    # Quality Metrics Section
    html += """
    <div class="section">
        <h2>Quality Metrics by MAF</h2>
        <div class="table-responsive">
            <table>
                <thead>
                    <tr>
                        <th>MAF Bin</th>
                        <th>Variant Count</th>
                        <th>Percentage</th>
                        <th>Mean R²</th>
                    </tr>
                </thead>
                <tbody>
"""
    
    # Add MAF distribution data
    maf_bins = genome_summary.get('maf_bins', {})
    r2_by_maf = genome_summary.get('mean_r2_by_maf', {})
    total_maf_variants = sum(maf_bins.values())
    
    for bin_name in sorted(maf_bins.keys()):
        count = maf_bins[bin_name]
        percentage = (count / total_maf_variants * 100) if total_maf_variants > 0 else 0
        mean_r2 = r2_by_maf.get(bin_name, {}).get('mean', 'N/A')
        r2_class = get_quality_class(mean_r2) if isinstance(mean_r2, float) else ''
        
        html += f"""
                    <tr>
                        <td>{bin_name}</td>
                        <td>{format_number(count)}</td>
                        <td>{percentage:.1f}%</td>
                        <td class="{r2_class}">{format_number(mean_r2) if mean_r2 != 'N/A' else 'N/A'}</td>
                    </tr>
"""
    
    html += """
                </tbody>
            </table>
        </div>
    </div>
"""
    
    # Chromosome Details Section
    html += """
    <div class="section">
        <h2>Chromosome-Level Summary</h2>
        <div class="chromosome-grid">
"""
    
    for chr_detail in genome_summary.get('chromosome_details', []):
        chr_name = chr_detail['chromosome']
        variants = chr_detail.get('variants', 0) or 0
        well_imputed = chr_detail.get('well_imputed', 0) or 0
        well_imputed_pct = (well_imputed / variants * 100) if variants > 0 else 0
        mean_r2 = chr_detail.get('mean_r2', 0) or 0
        
        html += f"""
            <div class="chromosome-card">
                <h4>Chromosome {chr_name}</h4>
                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                    <div>
                        <small>Chunks:</small> <strong>{chr_detail.get('chunks', 0)}</strong>
                    </div>
                    <div>
                        <small>Variants:</small> <strong>{format_number(variants)}</strong>
                    </div>
                    <div>
                        <small>Mean R²:</small> <strong class="{get_quality_class(mean_r2)}">{format_number(mean_r2)}</strong>
                    </div>
                    <div>
                        <small>Concordance:</small> <strong>{format_number(chr_detail.get('mean_concordance', 0) or 0)}</strong>
                    </div>
                </div>
                <div style="margin-top: 10px;">
                    <small>Well-imputed rate:</small>
                    <div class="progress-bar">
                        <div class="progress-fill" style="width: {well_imputed_pct:.1f}%"></div>
                    </div>
                    <small>{well_imputed_pct:.1f}% ({format_number(well_imputed)} / {format_number(variants)})</small>
                </div>
            </div>
"""
    
    html += """
        </div>
    </div>
"""
    
    # Plots Section (if available)
    if genome_plots_path and Path(genome_plots_path).exists():
        encoded_plots = encode_pdf_as_base64(genome_plots_path)
        if encoded_plots:
            html += f"""
    <div class="section">
        <h2>Visualization</h2>
        <div class="plot-container">
            <p>Comprehensive plots are available in the accompanying PDF file.</p>
            <embed src="data:application/pdf;base64,{encoded_plots}" width="100%" height="800px" type="application/pdf">
        </div>
    </div>
"""
    
    # Info Score Statistics (if available)
    if genome_summary.get('genome_info_stats'):
        info_stats = genome_summary['genome_info_stats']
        html += f"""
    <div class="section">
        <h2>Info Score Statistics</h2>
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="label">Mean Info Score</div>
                <div class="value">{format_number(info_stats.get('mean', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Median Info Score</div>
                <div class="value">{format_number(info_stats.get('median', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Min Info Score</div>
                <div class="value">{format_number(info_stats.get('min', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Max Info Score</div>
                <div class="value">{format_number(info_stats.get('max', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Variants Above 0.3</div>
                <div class="value">{format_number(info_stats.get('total_above_0.3', 0))}</div>
            </div>
            <div class="metric-card">
                <div class="label">Variants Above 0.8</div>
                <div class="value good">{format_number(info_stats.get('total_above_0.8', 0))}</div>
            </div>
        </div>
    </div>
"""
    
    # Footer
    html += f"""
    <div class="footer">
        <p>Generated by Genotype Imputation Pipeline v1.0</p>
        <p>Report generated on {timestamp}</p>
    </div>
</body>
</html>
"""
    
    return html


def get_quality_class(value):
    """Determine quality class for color coding."""
    if isinstance(value, (int, float)):
        if value >= 0.8:
            return 'good'
        elif value >= 0.3:
            return 'warning'
        else:
            return 'poor'
    return ''


def main():
    parser = argparse.ArgumentParser(description='Generate final imputation report')
    parser.add_argument('--genome-summary', required=True,
                       help='Genome summary JSON file')
    parser.add_argument('--genome-plots', required=True,
                       help='Genome plots PDF file')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    parser.add_argument('--population', default='Unknown',
                       help='Population')
    parser.add_argument('--study', default='Unknown',
                       help='Study name')
    
    args = parser.parse_args()
    
    # Load genome summary
    logger.info(f"Loading genome summary from {args.genome_summary}")
    genome_summary = load_genome_summary(args.genome_summary)
    
    # Generate HTML report
    logger.info("Generating HTML report...")
    html_content = generate_html_report(
        genome_summary,
        args.genome_plots,
        args.dataset,
        args.ref_name,
        args.population,
        args.study
    )
    
    # Save HTML report
    html_output = f"{args.output_prefix}_{args.ref_name}.final_report.html"
    with open(html_output, 'w') as f:
        f.write(html_content)
    logger.info(f"Saved HTML report to {html_output}")
    
    # Generate PDF report (if pdfkit is available)
    pdf_output = f"{args.output_prefix}_{args.ref_name}.final_report.pdf"
    
    pdf_created = False
    if PDFKIT_AVAILABLE:
        try:
            # Configure PDF options
            options = {
                'page-size': 'A4',
                'margin-top': '0.75in',
                'margin-right': '0.75in',
                'margin-bottom': '0.75in',
                'margin-left': '0.75in',
                'encoding': "UTF-8",
                'no-outline': None,
                'enable-local-file-access': None
            }
            
            # Try to generate PDF
            pdfkit.from_string(html_content, pdf_output, options=options)
            logger.info(f"Saved PDF report to {pdf_output}")
            pdf_created = True
        except Exception as e:
            logger.warning(f"Could not generate PDF report with pdfkit: {e}")
    
    if not pdf_created:
        logger.info("pdfkit not available, creating placeholder PDF with matplotlib...")
        # Create a simple placeholder PDF
        with PdfPages(pdf_output) as pdf:
            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.text(0.5, 0.5, f'Imputation Report\n\nDataset: {args.dataset}\nReference: {args.ref_name}\n\nPlease view the HTML report for full details.',
                   ha='center', va='center', fontsize=14, wrap=True)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
        logger.info(f"Created placeholder PDF at {pdf_output}")
    
    # Print summary
    print(f"\nFinal Report Generated:")
    print(f"  Dataset: {args.dataset}")
    print(f"  Reference: {args.ref_name}")
    print(f"  HTML Report: {html_output}")
    print(f"  PDF Report: {pdf_output}")
    print(f"  Total Variants: {genome_summary.get('total_variants', 0):,}")
    mean_r2_value = genome_summary.get('mean_r2', 0) or 0
    print(f"  Mean R²: {mean_r2_value:.4f}")
    
    logger.info("Report generation completed successfully")


if __name__ == '__main__':
    main()