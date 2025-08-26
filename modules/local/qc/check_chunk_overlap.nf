process CHECK_CHUNK_OVERLAP {
    tag "$chunk_name:$ref_name"
    label 'process_single'

    input:
    tuple val(meta), val(chunk_name), val(region), path(target_vcf), path(target_index), val(ref_name), path(ref_vcf)

    output:
    tuple val(meta), val(chunk_name), val(ref_name), path("${chunk_name}_${ref_name}_overlap.txt"), emit: overlap_report
    tuple val(meta), val(chunk_name), val(ref_name), path("${chunk_name}_${ref_name}_summary.txt"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${chunk_name}_${ref_name}"
    """
    echo "Checking overlap for chunk ${chunk_name} (region: ${region}) with reference panel ${ref_name}"
    
    # Ensure reference VCF has an index
    if [ ! -f "${ref_vcf}.tbi" ] && [ ! -f "${ref_vcf}.csi" ]; then
        echo "Creating index for reference panel..."
        bcftools index -t ${ref_vcf} --threads ${task.cpus}
    fi
    
    # Get variant positions from target chunk
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${target_vcf} | sort -k1,1 -k2,2n > target_variants.txt
    n_target=\$(wc -l < target_variants.txt)
    
    # Get variant positions from reference panel IN THE SAME REGION
    bcftools query -r ${region} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${ref_vcf} | sort -k1,1 -k2,2n > ref_variants.txt
    n_ref=\$(wc -l < ref_variants.txt)
    
    # Find exact matches (same position, ref, alt)
    comm -12 target_variants.txt ref_variants.txt > exact_matches.txt
    n_exact=\$(wc -l < exact_matches.txt)
    
    # Find position matches (regardless of alleles)
    cut -f1,2 target_variants.txt | sort -u > target_positions.txt
    cut -f1,2 ref_variants.txt | sort -u > ref_positions.txt
    comm -12 target_positions.txt ref_positions.txt > position_matches.txt
    n_pos_match=\$(wc -l < position_matches.txt)
    
    # Find target variants not in reference
    comm -23 target_variants.txt ref_variants.txt > target_only.txt
    n_target_only=\$(wc -l < target_only.txt)
    
    # Calculate overlap statistics
    if [ \$n_target -gt 0 ]; then
        overlap_pct=\$(echo "scale=2; 100 * \$n_exact / \$n_target" | bc)
        pos_overlap_pct=\$(echo "scale=2; 100 * \$n_pos_match / \$n_target" | bc)
    else
        overlap_pct=0
        pos_overlap_pct=0
    fi
    
    # Create detailed overlap report
    cat > ${prefix}_overlap.txt <<EOF
====================================
CHUNK OVERLAP REPORT
====================================
Chunk: ${chunk_name}
Reference Panel: ${ref_name}
------------------------------------

TARGET CHUNK STATISTICS:
  Total variants: \$n_target
  
REFERENCE CHUNK STATISTICS:
  Total variants: \$n_ref
  
OVERLAP STATISTICS:
  Exact matches (pos+ref+alt): \$n_exact (\${overlap_pct}%)
  Position matches: \$n_pos_match (\${pos_overlap_pct}%)
  Target-only variants: \$n_target_only
  
OVERLAP QUALITY:
EOF
    
    # Assess overlap quality
    if [ \$(echo "\$overlap_pct >= 80" | bc -l) -eq 1 ]; then
        echo "  Status: EXCELLENT (≥80% exact match)" >> ${prefix}_overlap.txt
    elif [ \$(echo "\$overlap_pct >= 60" | bc -l) -eq 1 ]; then
        echo "  Status: GOOD (≥60% exact match)" >> ${prefix}_overlap.txt
    elif [ \$(echo "\$overlap_pct >= 40" | bc -l) -eq 1 ]; then
        echo "  Status: ACCEPTABLE (≥40% exact match)" >> ${prefix}_overlap.txt
    elif [ \$(echo "\$overlap_pct >= 20" | bc -l) -eq 1 ]; then
        echo "  Status: POOR (≥20% exact match)" >> ${prefix}_overlap.txt
    else
        echo "  Status: VERY POOR (<20% exact match)" >> ${prefix}_overlap.txt
    fi
    
    # Add recommendations
    echo "" >> ${prefix}_overlap.txt
    echo "RECOMMENDATIONS:" >> ${prefix}_overlap.txt
    if [ \$(echo "\$overlap_pct < 40" | bc -l) -eq 1 ]; then
        echo "  ⚠ Low overlap detected. Consider:" >> ${prefix}_overlap.txt
        echo "    - Using a different reference panel" >> ${prefix}_overlap.txt
        echo "    - Checking genome build compatibility" >> ${prefix}_overlap.txt
        echo "    - Verifying chromosome naming conventions" >> ${prefix}_overlap.txt
    else
        echo "  ✓ Overlap is sufficient for imputation" >> ${prefix}_overlap.txt
    fi
    
    # Sample of non-matching variants for debugging
    if [ \$n_target_only -gt 0 ]; then
        echo "" >> ${prefix}_overlap.txt
        echo "SAMPLE OF TARGET-ONLY VARIANTS (first 10):" >> ${prefix}_overlap.txt
        head -10 target_only.txt | while read chr pos ref alt; do
            echo "  \$chr:\$pos \$ref>\$alt" >> ${prefix}_overlap.txt
        done
        if [ \$n_target_only -gt 10 ]; then
            echo "  ... and \$((n_target_only - 10)) more" >> ${prefix}_overlap.txt
        fi
    fi
    
    echo "====================================
    " >> ${prefix}_overlap.txt
    
    # Create summary file for aggregation
    echo "${chunk_name}\\t${ref_name}\\t\$n_target\\t\$n_ref\\t\$n_exact\\t\${overlap_pct}\\t\$n_pos_match\\t\${pos_overlap_pct}" > ${prefix}_summary.txt
    
    # Display summary
    echo "Overlap: \$n_exact/\$n_target variants (\${overlap_pct}%) match exactly with reference"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}