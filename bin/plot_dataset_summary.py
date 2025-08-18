#!/usr/bin/env python3
"""
Create summary visualization for dataset-level metrics
"""

import argparse
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path


def create_summary_plot(summary_data, output_file):
    """Create a comprehensive summary plot for dataset-level metrics"""
    
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f"Dataset Summary: {summary_data['dataset_id']}\n"
                 f"Reference Panel: {summary_data['reference_panel']}", 
                 fontsize=16, fontweight='bold')
    
    # Create grid for subplots
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Chunks processed (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(0.5, 0.5, f"{summary_data['chunks_processed']:,}", 
             ha='center', va='center', fontsize=24, fontweight='bold')
    ax1.text(0.5, 0.2, "Chunks Processed", 
             ha='center', va='center', fontsize=12)
    if summary_data.get('chunks_failed', 0) > 0:
        ax1.text(0.5, 0.05, f"({summary_data['chunks_failed']} failed)", 
                 ha='center', va='center', fontsize=10, color='red')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis('off')
    
    # 2. Variant counts (top middle)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.text(0.5, 0.7, f"{summary_data['total_variants']:,}", 
             ha='center', va='center', fontsize=20, fontweight='bold')
    ax2.text(0.5, 0.5, "Total Variants", 
             ha='center', va='center', fontsize=12)
    ax2.text(0.5, 0.3, f"{summary_data['well_imputed_variants']:,}", 
             ha='center', va='center', fontsize=16, color='green')
    ax2.text(0.5, 0.15, "Well Imputed", 
             ha='center', va='center', fontsize=10, color='green')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis('off')
    
    # 3. Imputation rate gauge (top right)
    ax3 = fig.add_subplot(gs[0, 2])
    imputation_rate = summary_data.get('imputation_rate', 0)
    
    # Create gauge chart
    theta = np.linspace(0, np.pi, 100)
    r_inner = 0.7
    r_outer = 1.0
    
    # Background arc
    ax3.fill_between(theta, r_inner, r_outer, transform=ax3.transData, 
                     color='lightgray', alpha=0.5)
    
    # Filled arc based on imputation rate
    theta_fill = np.linspace(0, np.pi * imputation_rate / 100, 100)
    color = 'green' if imputation_rate > 80 else 'orange' if imputation_rate > 60 else 'red'
    ax3.fill_between(theta_fill, r_inner, r_outer, transform=ax3.transData, 
                     color=color, alpha=0.8)
    
    ax3.text(0, -0.2, f"{imputation_rate:.1f}%", 
             ha='center', va='center', fontsize=20, fontweight='bold')
    ax3.text(0, -0.4, "Imputation Rate", 
             ha='center', va='center', fontsize=12)
    ax3.set_xlim(-1.2, 1.2)
    ax3.set_ylim(-0.5, 1.2)
    ax3.axis('off')
    
    # 4. R² distribution (middle row, full width)
    ax4 = fig.add_subplot(gs[1, :])
    
    # Create bar chart for R² statistics
    r2_stats = {
        'Mean': summary_data.get('mean_r2', 0),
        'Median': summary_data.get('median_r2', 0) if summary_data.get('median_r2') is not None else summary_data.get('mean_r2', 0),
        '25th %ile': summary_data.get('percentile_25_r2', 0) if summary_data.get('percentile_25_r2') is not None else summary_data.get('mean_r2', 0),
        '75th %ile': summary_data.get('percentile_75_r2', 0) if summary_data.get('percentile_75_r2') is not None else summary_data.get('mean_r2', 0)
    }
    
    positions = np.arange(len(r2_stats))
    values = list(r2_stats.values())
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
    
    bars = ax4.bar(positions, values, color=colors, alpha=0.8)
    ax4.set_xticks(positions)
    ax4.set_xticklabels(r2_stats.keys())
    ax4.set_ylabel('R² Value', fontsize=12)
    ax4.set_title('R² Statistics Across All Chunks', fontsize=12, fontweight='bold')
    ax4.set_ylim(0, 1)
    ax4.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10)
    
    # Add reference line at 0.8 (typical good R² threshold)
    ax4.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='Good R² (0.8)')
    ax4.legend(loc='upper right')
    
    # 5. MAF summary (bottom left)
    ax5 = fig.add_subplot(gs[2, 0])
    mean_maf = summary_data.get('mean_maf', 0)
    ax5.text(0.5, 0.5, f"{mean_maf:.4f}", 
             ha='center', va='center', fontsize=20, fontweight='bold')
    ax5.text(0.5, 0.2, "Mean MAF", 
             ha='center', va='center', fontsize=12)
    ax5.set_xlim(0, 1)
    ax5.set_ylim(0, 1)
    ax5.axis('off')
    
    # 6. Quality summary (bottom middle and right)
    ax6 = fig.add_subplot(gs[2, 1:])
    
    # Create quality assessment text
    quality_text = []
    
    # Assess overall quality
    if summary_data.get('mean_r2', 0) > 0.8:
        quality_text.append(("• Excellent imputation quality (R² > 0.8)", 'green'))
    elif summary_data.get('mean_r2', 0) > 0.6:
        quality_text.append(("• Good imputation quality (R² > 0.6)", 'orange'))
    else:
        quality_text.append(("• Low imputation quality (R² < 0.6)", 'red'))
    
    if imputation_rate > 90:
        quality_text.append(("• High variant coverage (>90%)", 'green'))
    elif imputation_rate > 70:
        quality_text.append(("• Moderate variant coverage (70-90%)", 'orange'))
    else:
        quality_text.append(("• Low variant coverage (<70%)", 'red'))
    
    if summary_data.get('chunks_failed', 0) == 0:
        quality_text.append(("• All chunks processed successfully", 'green'))
    else:
        quality_text.append((f"• {summary_data['chunks_failed']} chunks failed", 'red'))
    
    # Display quality assessment
    y_pos = 0.8
    for text, color in quality_text:
        ax6.text(0.05, y_pos, text, fontsize=11, color=color, 
                transform=ax6.transAxes, fontweight='bold')
        y_pos -= 0.25
    
    ax6.set_title('Quality Assessment', fontsize=12, fontweight='bold')
    ax6.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Dataset summary plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Create dataset-level summary visualization')
    parser.add_argument('--summary-json', required=True, 
                       help='JSON file with aggregated summary statistics')
    parser.add_argument('--output-prefix', required=True,
                       help='Output prefix for plot files')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    
    args = parser.parse_args()
    
    # Load summary data
    with open(args.summary_json, 'r') as f:
        summary_data = json.load(f)
    
    # Create output filename
    output_file = f"{args.output_prefix}_{args.ref_name}.dataset_summary.png"
    
    # Create summary plot
    create_summary_plot(summary_data, output_file)


if __name__ == '__main__':
    main()