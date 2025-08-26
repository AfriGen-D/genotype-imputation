process GENERATE_FINAL_REPORT {
    tag "$meta.id"
    label 'process_low'
    label 'latex_reporting'
    container 'mamana/chipimputation-latex-reporting:1.0.0'
    
    publishDir "${params.outdir}/${meta.dataset ?: meta.id}/${ref_name}/reports", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(genome_plots, stageAs: 'plots/genome/*'), path(chr_plots, stageAs: 'plots/chr/*'), path(chunk_plots, stageAs: 'plots/chunk/*'), path(stats_files, stageAs: 'stats/*')
    
    output:
    tuple val(meta), val(ref_name), path("*.final_report.pdf"), emit: report
    tuple val(meta), val(ref_name), path("*.final_report.tex"), emit: latex
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dataset = meta.dataset ?: meta.id
    def report_base = "${dataset}_${ref_name}.final_report"
    """
    # Create directory structure
    mkdir -p output_dir/{plots,stats}
    mkdir -p output_dir/plots/{genome_level,chromosome_level,chunk_level}
    mkdir -p output_dir/stats/{genome_stats,chromosome_stats,chunk_stats}
    
    # Organize plots
    if [ -d plots/genome ]; then
        cp -r plots/genome/* output_dir/plots/genome_level/ 2>/dev/null || true
    fi
    if [ -d plots/chr ]; then
        cp -r plots/chr/* output_dir/plots/chromosome_level/ 2>/dev/null || true
    fi
    if [ -d plots/chunk ]; then
        cp -r plots/chunk/* output_dir/plots/chunk_level/ 2>/dev/null || true
    fi
    
    # Organize statistics
    if [ -d stats ]; then
        find stats -name "*genome*.json" -exec cp {} output_dir/stats/genome_stats/ \\; 2>/dev/null || true
        find stats -name "*chr*.json" -exec cp {} output_dir/stats/chromosome_stats/ \\; 2>/dev/null || true
        find stats -name "*chunk*.json" -exec cp {} output_dir/stats/chunk_stats/ \\; 2>/dev/null || true
    fi
    
    # Copy report generator scripts
    cp ${projectDir}/bin/generate_latex_report.py .
    cp ${projectDir}/bin/generate_pdf_report.py .
    
    # Generate LaTeX report
    echo "Generating LaTeX report..."
    python3 generate_latex_report.py \\
        --dataset "${dataset}" \\
        --ref-panel "${ref_name}" \\
        --output-dir output_dir \\
        --output ${report_base}.tex
    
    # Compile LaTeX to PDF (should work in latex container)
    echo "Compiling LaTeX to PDF..."
    if command -v pdflatex &> /dev/null; then
        # Use the compile-report helper if available
        if command -v compile-report &> /dev/null; then
            compile-report ${report_base}.tex
        else
            # Manual compilation
            pdflatex -interaction=nonstopmode ${report_base}.tex || true
            pdflatex -interaction=nonstopmode ${report_base}.tex || true
            
            # Clean up auxiliary files
            rm -f ${report_base}.{aux,log,out,toc} 2>/dev/null || true
        fi
    fi
    
    # If LaTeX PDF failed, use Python PDF generator as fallback
    if [ ! -f ${report_base}.pdf ]; then
        echo "LaTeX compilation failed, using Python PDF generator..."
        python3 generate_pdf_report.py \\
            --dataset "${dataset}" \\
            --ref-panel "${ref_name}" \\
            --output-dir output_dir \\
            --output ${report_base}.pdf
    fi
    
    # Verify PDF was created
    if [ -f ${report_base}.pdf ]; then
        echo "PDF report successfully generated: ${report_base}.pdf"
        ls -lh ${report_base}.pdf
    else
        echo "Warning: PDF generation failed, but LaTeX source is available"
    fi
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "N/A")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)" 2>/dev/null || echo "N/A")
        latex: \$(pdflatex --version 2>/dev/null | head -1 || echo "Not installed")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dataset = meta.dataset ?: meta.id
    def report_base = "${dataset}_${ref_name}.final_report"
    """
    touch ${report_base}.pdf
    touch ${report_base}.tex
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
        matplotlib: 3.7.0
        latex: pdfTeX 3.14159265
    END_VERSIONS
    """
}