process AGGREGATE_OVERLAP_REPORTS {
    tag "overlap_summary"
    label 'process_single'
    publishDir "${params.outdir}/qc/overlap_reports", mode: params.publish_dir_mode
    
    input:
    path overlap_reports
    path summary_files

    output:
    path "chunk_overlap_summary.txt", emit: summary
    path "chunk_overlap_detailed.txt", emit: detailed
    path "chunk_overlap_matrix.csv", emit: matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Aggregating overlap reports from all chunks..."
    
    # Create header for summary
    cat > chunk_overlap_summary.txt <<EOF
====================================
CHUNK OVERLAP SUMMARY REPORT
====================================
Generated: \$(date)
Total chunks analyzed: \$(ls *_summary.txt 2>/dev/null | wc -l)
====================================

EOF
    
    # Create CSV header for matrix
    echo "Chunk,Reference_Panel,Target_Variants,Ref_Variants,Exact_Matches,Exact_Match_%,Position_Matches,Position_Match_%" > chunk_overlap_matrix.csv
    
    # Process summary files
    for summary in *_summary.txt; do
        cat \$summary >> chunk_overlap_matrix.csv
    done
    
    # Calculate overall statistics
    total_chunks=0
    excellent=0
    good=0
    acceptable=0
    poor=0
    very_poor=0
    
    # Analyze each overlap report
    for report in *_overlap.txt; do
        if [ -f "\$report" ]; then
            total_chunks=\$((total_chunks + 1))
            
            # Extract status from report
            status=\$(grep "Status:" \$report | head -1 | cut -d: -f2 | cut -d'(' -f1 | xargs)
            
            case "\$status" in
                "EXCELLENT") excellent=\$((excellent + 1)) ;;
                "GOOD") good=\$((good + 1)) ;;
                "ACCEPTABLE") acceptable=\$((acceptable + 1)) ;;
                "POOR") poor=\$((poor + 1)) ;;
                "VERY POOR") very_poor=\$((very_poor + 1)) ;;
            esac
        fi
    done
    
    # Add statistics to summary
    cat >> chunk_overlap_summary.txt <<EOF
OVERALL STATISTICS:
  Excellent (≥80%): \$excellent chunks
  Good (≥60%): \$good chunks
  Acceptable (≥40%): \$acceptable chunks
  Poor (≥20%): \$poor chunks
  Very Poor (<20%): \$very_poor chunks

EOF
    
    # Calculate average overlap across all chunks
    if [ -f chunk_overlap_matrix.csv ]; then
        avg_overlap=\$(tail -n +2 chunk_overlap_matrix.csv | awk -F, '{sum+=\$6; count++} END {if(count>0) print sum/count; else print 0}')
        echo "Average exact match overlap: \${avg_overlap}%" >> chunk_overlap_summary.txt
    fi
    
    # Identify problematic chunks
    echo "" >> chunk_overlap_summary.txt
    echo "CHUNKS REQUIRING ATTENTION (overlap <40%):" >> chunk_overlap_summary.txt
    tail -n +2 chunk_overlap_matrix.csv | awk -F, '\$6 < 40 {print "  - " \$1 " with " \$2 ": " \$6 "%"}' >> chunk_overlap_summary.txt
    
    # Combine all detailed reports
    echo "====================================
DETAILED CHUNK OVERLAP REPORTS
====================================
" > chunk_overlap_detailed.txt
    
    for report in *_overlap.txt; do
        if [ -f "\$report" ]; then
            cat \$report >> chunk_overlap_detailed.txt
            echo "" >> chunk_overlap_detailed.txt
        fi
    done
    
    # Display summary
    echo "Overlap analysis complete:"
    echo "  Excellent: \$excellent chunks"
    echo "  Good: \$good chunks"
    echo "  Acceptable: \$acceptable chunks"
    echo "  Poor: \$poor chunks"
    echo "  Very Poor: \$very_poor chunks"
    
    if [ \$poor -gt 0 ] || [ \$very_poor -gt 0 ]; then
        echo ""
        echo "⚠ WARNING: \$((poor + very_poor)) chunks have low overlap (<40%)"
        echo "Check chunk_overlap_summary.txt for details"
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}