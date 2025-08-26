#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate comprehensive LaTeX-based imputation report with plots and interpretations
"""

import os
import sys
import json
import glob
import argparse
import subprocess
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np

def escape_latex(text):
    """Escape special LaTeX characters"""
    replacements = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text

def load_genome_statistics(stats_dir):
    """Load genome-wide statistics"""
    stats = {
        'total_variants': 0,
        'mean_r2': 0.0,
        'mean_info': 0.0,
        'well_imputed_pct': 0.0,
        'chromosomes_processed': 0,
        'total_chunks': 0
    }
    
    genome_stats_files = glob.glob(f"{stats_dir}/genome_stats/*.summary.json")
    if genome_stats_files:
        with open(genome_stats_files[0]) as f:
            loaded_stats = json.load(f)
            stats.update(loaded_stats)
    
    # Count chunks
    chunk_files = glob.glob(f"{stats_dir}/chunk_stats/*.json")
    stats['total_chunks'] = len(chunk_files)
    
    return stats

def load_chromosome_statistics(stats_dir):
    """Load chromosome-level statistics"""
    chr_stats = []
    for chr_num in range(1, 23):
        chr_files = glob.glob(f"{stats_dir}/chromosome_stats/*chr{chr_num}*.summary.json")
        if chr_files:
            with open(chr_files[0]) as f:
                stats = json.load(f)
                stats['chromosome'] = chr_num
                chr_stats.append(stats)
    return chr_stats

def interpret_results(stats):
    """Generate interpretation of results"""
    interpretations = []
    
    mean_r2 = stats.get('mean_r2', 0)
    mean_info = stats.get('mean_info', 0)
    well_imputed = stats.get('well_imputed_pct', 0)
    
    # Overall quality assessment
    if mean_r2 >= 0.8:
        quality = "excellent"
        interpretations.append("The overall imputation quality is excellent with high accuracy across the genome.")
    elif mean_r2 >= 0.6:
        quality = "good"
        interpretations.append("The imputation shows good overall quality suitable for most downstream analyses.")
    elif mean_r2 >= 0.4:
        quality = "moderate"
        interpretations.append("The imputation quality is moderate. Consider additional QC filtering for association studies.")
    else:
        quality = "poor"
        interpretations.append("The imputation quality is below optimal. Review input data quality and reference panel compatibility.")
    
    # INFO score interpretation
    if mean_info >= 0.8:
        interpretations.append(f"The average INFO score of {mean_info:.3f} indicates high confidence in imputed genotypes.")
    elif mean_info >= 0.5:
        interpretations.append(f"The average INFO score of {mean_info:.3f} suggests moderate imputation certainty.")
    else:
        interpretations.append(f"The low average INFO score of {mean_info:.3f} suggests uncertainty in imputed genotypes.")
    
    # Well-imputed variants
    if well_imputed >= 80:
        interpretations.append(f"Approximately {well_imputed:.1f}\\% of variants are well-imputed (R² ≥ 0.3), indicating good coverage.")
    elif well_imputed >= 60:
        interpretations.append(f"About {well_imputed:.1f}\\% of variants meet the quality threshold (R² ≥ 0.3).")
    else:
        interpretations.append(f"Only {well_imputed:.1f}\\% of variants are well-imputed, suggesting limited imputation effectiveness.")
    
    return interpretations

def generate_latex_report(dataset_id, ref_panel, output_dir, output_tex):
    """Generate LaTeX report"""
    
    stats_dir = f"{output_dir}/stats"
    plots_dir = f"{output_dir}/plots"
    
    # Load statistics
    genome_stats = load_genome_statistics(stats_dir)
    chr_stats = load_chromosome_statistics(stats_dir)
    interpretations = interpret_results(genome_stats)
    
    # Start LaTeX document
    latex_content = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry}
\usepackage{graphicx}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{xcolor}
\usepackage{colortbl}
\usepackage{hyperref}
\usepackage{fancyhdr}
\usepackage{lastpage}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=1.17}

\geometry{
    left=25mm,
    right=25mm,
    top=30mm,
    bottom=30mm
}

\definecolor{chipblue}{RGB}{44, 62, 80}
\definecolor{chipgreen}{RGB}{39, 174, 96}
\definecolor{chiporange}{RGB}{230, 126, 34}
\definecolor{chipred}{RGB}{192, 57, 43}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small ChiPImputation Report}
\fancyhead[R]{\small \today}
\fancyfoot[C]{\small Page \thepage\ of \pageref{LastPage}}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}

\title{
    \vspace{-2cm}
    \huge \textbf{ChiPImputation}\\
    \vspace{0.5cm}
    \Large Comprehensive Genotype Imputation Report\\
    \vspace{0.3cm}
    \large Dataset: """ + escape_latex(dataset_id) + r"""\\
    Reference Panel: """ + escape_latex(ref_panel) + r"""
}
\author{Generated by ChiPImputation Pipeline v1.0.0}
\date{\today}

\begin{document}

\maketitle
\thispagestyle{empty}
\vspace{2cm}

\begin{abstract}
This report presents comprehensive results from genotype imputation performed using the ChiPImputation pipeline. 
The analysis processed """ + f"{genome_stats['total_chunks']}" + r""" genomic chunks across """ + f"{genome_stats.get('chromosomes_processed', 22)}" + r""" chromosomes, 
imputing a total of """ + f"{genome_stats['total_variants']:,}" + r""" variants. 
The overall imputation quality achieved a mean R$^2$ of """ + f"{genome_stats['mean_r2']:.3f}" + r""" with """ + f"{genome_stats['well_imputed_pct']:.1f}" + r"""\% of variants 
meeting the quality threshold (R$^2$ $\geq$ 0.3).
\end{abstract}

\newpage
\tableofcontents
\newpage

\section{Executive Summary}

\subsection{Key Metrics}
\begin{table}[H]
\centering
\caption{Overall Imputation Performance Metrics}
\begin{tabular}{lc}
\toprule
\textbf{Metric} & \textbf{Value} \\
\midrule
Total Variants Imputed & """ + f"{genome_stats['total_variants']:,}" + r""" \\
Mean Imputation R$^2$ & """ + f"{genome_stats['mean_r2']:.3f}" + r""" \\
Mean INFO Score & """ + f"{genome_stats['mean_info']:.3f}" + r""" \\
Well-Imputed Variants (R$^2$ $\geq$ 0.3) & """ + f"{genome_stats['well_imputed_pct']:.1f}" + r"""\% \\
Chromosomes Processed & """ + f"{genome_stats.get('chromosomes_processed', 22)}" + r""" \\
Total Genomic Chunks & """ + f"{genome_stats['total_chunks']}" + r""" \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Results Interpretation}
\begin{itemize}
"""
    
    # Add interpretations
    for interpretation in interpretations:
        latex_content += f"    \\item {interpretation}\n"
    
    latex_content += r"""
\end{itemize}

\section{Genome-Wide Analysis}

\subsection{Imputation Accuracy Overview}
"""
    
    # Add genome-wide accuracy plot if it exists
    accuracy_plots = glob.glob(f"{plots_dir}/genome_level/*/genome_accuracy.pdf")
    if accuracy_plots:
        latex_content += r"""
\begin{figure}[H]
    \centering
    \includegraphics[width=\textwidth]{""" + accuracy_plots[0] + r"""}
    \caption{Genome-wide imputation accuracy metrics showing R$^2$, INFO scores, and well-imputed percentages across all chromosomes.}
    \label{fig:genome_accuracy}
\end{figure}

Figure~\ref{fig:genome_accuracy} presents the comprehensive accuracy metrics across the genome. 
The bar charts show chromosome-specific performance, while the scatter plot reveals the correlation 
between INFO scores and R$^2$ values, indicating consistency in imputation quality measures.
"""

    latex_content += r"""
\subsection{MAF vs R$^2$ Distribution}
"""
    
    # Add MAF vs R2 plot if it exists
    maf_plots = glob.glob(f"{plots_dir}/genome_level/*/genome_maf_r2.pdf")
    if maf_plots:
        latex_content += r"""
\begin{figure}[H]
    \centering
    \includegraphics[width=0.8\textwidth]{""" + maf_plots[0] + r"""}
    \caption{Relationship between minor allele frequency (MAF) and imputation accuracy (R$^2$).}
    \label{fig:maf_r2}
\end{figure}

The relationship between MAF and imputation accuracy (Figure~\ref{fig:maf_r2}) shows the expected pattern 
where common variants (higher MAF) typically achieve better imputation accuracy than rare variants. 
This is particularly important for association studies focusing on different frequency spectra.
"""

    latex_content += r"""
\subsection{Chromosome-Level Summary}

\begin{table}[H]
\centering
\caption{Imputation Performance by Chromosome}
\small
\begin{tabular}{cccccc}
\toprule
\textbf{Chr} & \textbf{Variants} & \textbf{Mean R$^2$} & \textbf{Mean INFO} & \textbf{Well-Imputed (\%)} \\
\midrule
"""
    
    # Add chromosome statistics
    for chr_stat in sorted(chr_stats, key=lambda x: x['chromosome']):
        latex_content += f"{chr_stat['chromosome']} & "
        latex_content += f"{chr_stat.get('n_variants', 0):,} & "
        latex_content += f"{chr_stat.get('mean_r2', 0):.3f} & "
        latex_content += f"{chr_stat.get('mean_info', 0):.3f} & "
        latex_content += f"{chr_stat.get('well_imputed_pct', 0):.1f} \\\\\n"
    
    latex_content += r"""
\bottomrule
\end{tabular}
\end{table}

\section{Quality Control}

\subsection{Pre-Imputation QC}
The following quality control steps were applied before imputation:
\begin{itemize}
    \item \textbf{Duplicate Removal:} Identification and removal of duplicate variants
    \item \textbf{Multi-allelic Splitting:} Splitting of multi-allelic variants into biallelic records
    \item \textbf{Minor Allele Count:} Filtering variants with MAC $\geq$ 1
    \item \textbf{Site Missingness:} Maximum allowed missingness of 5\%
    \item \textbf{Hardy-Weinberg Equilibrium:} HWE p-value threshold of $10^{-5}$
\end{itemize}

\subsection{Post-Imputation QC}
Post-imputation quality metrics and thresholds:
\begin{itemize}
    \item \textbf{INFO Score Threshold:} $\geq$ 0.3 for inclusion in analysis
    \item \textbf{R$^2$ Threshold:} $\geq$ 0.3 for well-imputed classification
    \item \textbf{MAF Threshold:} $\geq$ 0.01 for downstream analyses
\end{itemize}

\section{Methodology}

\subsection{Pipeline Configuration}
\begin{table}[H]
\centering
\caption{Imputation Pipeline Parameters}
\begin{tabular}{ll}
\toprule
\textbf{Parameter} & \textbf{Value} \\
\midrule
Reference Panel & """ + escape_latex(ref_panel) + r""" \\
Genome Build & hg38/GRCh38 \\
Phasing Algorithm & EAGLE v2.4.1 \\
Imputation Algorithm & MINIMAC4 v4.1.2 \\
Chunk Size & 25 Mb \\
Buffer Size & 1 Mb \\
Effective Population Size & 20,000 \\
MCMC Iterations & 10 \\
Burn-in Iterations & 2 \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Computational Resources}
The imputation was performed using high-performance computing infrastructure with:
\begin{itemize}
    \item SLURM job scheduler for parallel processing
    \item Singularity containers for reproducibility
    \item Maximum 50GB memory per process
    \item Up to 500 parallel jobs
\end{itemize}

\section{Recommendations}

Based on the imputation results, we recommend the following for downstream analyses:

\begin{enumerate}
    \item \textbf{Variant Filtering:} Apply INFO score threshold of $\geq$ 0.8 for association studies
    \item \textbf{MAF-Specific Thresholds:} Consider using MAF-specific R$^2$ thresholds:
        \begin{itemize}
            \item MAF $\geq$ 0.05: R$^2$ $\geq$ 0.3
            \item 0.01 $\leq$ MAF < 0.05: R$^2$ $\geq$ 0.6
            \item MAF < 0.01: R$^2$ $\geq$ 0.8
        \end{itemize}
    \item \textbf{Quality Review:} Chromosomes with mean R$^2$ < 0.5 should be reviewed for data quality issues
    \item \textbf{Population Stratification:} Consider population-specific imputation if mixed ancestry is detected
\end{enumerate}

\section{Additional Plots and Visualizations}
"""
    
    # Add chromosome-level plots section
    chr_perf_plots = glob.glob(f"{plots_dir}/chromosome_level/chr_performance/*.pdf")[:3]  # Show first 3
    if chr_perf_plots:
        latex_content += r"""
\subsection{Chromosome-Specific Performance}
Selected chromosome-level performance visualizations:
"""
        for i, plot in enumerate(chr_perf_plots, 1):
            chr_num = Path(plot).stem.split('_')[2].replace('chr', '')
            latex_content += r"""
\begin{figure}[H]
    \centering
    \includegraphics[width=0.9\textwidth]{""" + plot + r"""}
    \caption{Detailed performance metrics for chromosome """ + chr_num + r""".}
\end{figure}
"""
    
    latex_content += r"""
\section{Appendix}

\subsection{File Locations}
\begin{itemize}
    \item \textbf{Imputed VCFs:} \texttt{""" + escape_latex(f"{output_dir}/imputed/") + r"""}
    \item \textbf{Quality Metrics:} \texttt{""" + escape_latex(f"{output_dir}/stats/") + r"""}
    \item \textbf{Plots:} \texttt{""" + escape_latex(f"{output_dir}/plots/") + r"""}
    \item \textbf{Logs:} \texttt{""" + escape_latex(f"{output_dir}/logs/") + r"""}
\end{itemize}

\subsection{References}
\begin{enumerate}
    \item Loh, P.R., et al. (2016). Reference-based phasing using the Haplotype Reference Consortium panel. \textit{Nature Genetics}, 48(11), 1443-1448.
    \item Das, S., et al. (2016). Next-generation genotype imputation service and methods. \textit{Nature Genetics}, 48(10), 1284-1287.
    \item H3Africa Consortium (2014). Enabling the genomic revolution in Africa. \textit{Science}, 344(6190), 1346-1348.
\end{enumerate}

\end{document}
"""
    
    # Write LaTeX file
    with open(output_tex, 'w') as f:
        f.write(latex_content)
    
    print(f"LaTeX report generated: {output_tex}")
    
    # Try to compile to PDF if pdflatex is available
    try:
        pdf_file = output_tex.replace('.tex', '.pdf')
        subprocess.run(['pdflatex', '-interaction=nonstopmode', output_tex], 
                      capture_output=True, cwd=os.path.dirname(output_tex))
        # Run twice for references
        subprocess.run(['pdflatex', '-interaction=nonstopmode', output_tex], 
                      capture_output=True, cwd=os.path.dirname(output_tex))
        print(f"PDF report generated: {pdf_file}")
        
        # Clean up auxiliary files
        for ext in ['.aux', '.log', '.out', '.toc']:
            aux_file = output_tex.replace('.tex', ext)
            if os.path.exists(aux_file):
                os.remove(aux_file)
    except Exception as e:
        print(f"Could not compile PDF (pdflatex may not be installed): {e}")
        print("LaTeX source file is available for manual compilation")

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive LaTeX imputation report')
    parser.add_argument('--dataset', required=True, help='Dataset ID')
    parser.add_argument('--ref-panel', required=True, help='Reference panel name')
    parser.add_argument('--output-dir', required=True, help='Output directory with results')
    parser.add_argument('--output', required=True, help='Output LaTeX/PDF file')
    
    args = parser.parse_args()
    
    # Ensure output ends with .tex
    if not args.output.endswith('.tex'):
        args.output = args.output.replace('.pdf', '.tex')
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)
    
    generate_latex_report(
        args.dataset,
        args.ref_panel,
        args.output_dir,
        args.output
    )

if __name__ == '__main__':
    main()