process CHECK_HWE {
    tag "$meta.id"
    label 'process_medium'
    
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    val hwe_p_threshold
    
    output:
    tuple val(meta), path("${prefix}.hwe_stats.txt"), emit: stats
    tuple val(meta), path("${prefix}.hwe_violations.txt"), emit: violations  
    tuple val(meta), path("${prefix}.hwe_summary.txt"), emit: report
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pval_threshold = hwe_p_threshold ?: "0.001"
    """
    # Calculate Hardy-Weinberg equilibrium statistics
    bcftools +fill-tags ${vcf} -- -t HWE | \\
        bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%HWE\\n' > temp_hwe.txt
    
    # Process HWE results
    awk -v threshold=${pval_threshold} '
    BEGIN { 
        OFS="\\t"
        print "Chromosome\\tPosition\\tVariant_ID\\tREF\\tALT\\tHWE_P_Value\\tSignificance" > "${prefix}.hwe_stats.txt"
        print "Chromosome\\tPosition\\tVariant_ID\\tREF\\tALT\\tHWE_P_Value" > "${prefix}.hwe_violations.txt"
    }
    {
        variant_id = (\$3 != "." && \$3 != "") ? \$3 : \$1":"\$2
        hwe_p = (\$6 != "." && \$6 != "") ? \$6 : "NA"
        
        if(hwe_p != "NA") {
            significance = (hwe_p+0 < threshold+0) ? "VIOLATION" : "PASS"
            print \$1, \$2, variant_id, \$4, \$5, hwe_p, significance >> "${prefix}.hwe_stats.txt"
            
            if(hwe_p+0 < threshold+0) {
                print \$1, \$2, variant_id, \$4, \$5, hwe_p >> "${prefix}.hwe_violations.txt"
            }
        }
    }' temp_hwe.txt
    
    # Generate summary report
    echo "Hardy-Weinberg Equilibrium Analysis for ${meta.id}" > ${prefix}.hwe_summary.txt
    echo "Generated on: \$(date)" >> ${prefix}.hwe_summary.txt
    echo "P-value threshold: ${pval_threshold}" >> ${prefix}.hwe_summary.txt
    echo "" >> ${prefix}.hwe_summary.txt
    
    # Calculate summary statistics
    total_variants=\$(tail -n +2 ${prefix}.hwe_stats.txt | wc -l)
    violations=\$(tail -n +2 ${prefix}.hwe_violations.txt | wc -l)
    
    echo "SUMMARY STATISTICS:" >> ${prefix}.hwe_summary.txt
    echo "Total variants tested: \$total_variants" >> ${prefix}.hwe_summary.txt
    echo "HWE violations (p < ${pval_threshold}): \$violations" >> ${prefix}.hwe_summary.txt
    
    if [ \$total_variants -gt 0 ]; then
        violation_rate=\$(echo "scale=4; 100 * \$violations / \$total_variants" | bc -l)
        echo "Violation rate: \$violation_rate%" >> ${prefix}.hwe_summary.txt
    fi
    
    echo "" >> ${prefix}.hwe_summary.txt
    
    # P-value distribution
    if [ \$total_variants -gt 0 ]; then
        echo "P-VALUE DISTRIBUTION:" >> ${prefix}.hwe_summary.txt
        tail -n +2 ${prefix}.hwe_stats.txt | awk -F'\\t' '
        \$6 != "NA" {
            pval = \$6 + 0
            if(pval < 0.001) bin1++
            else if(pval < 0.01) bin2++  
            else if(pval < 0.05) bin3++
            else if(pval < 0.1) bin4++
            else bin5++
            total++
        }
        END {
            if(total > 0) {
                printf "p < 0.001: %d (%.2f%%)\\n", bin1+0, 100*(bin1+0)/total
                printf "0.001 ≤ p < 0.01: %d (%.2f%%)\\n", bin2+0, 100*(bin2+0)/total  
                printf "0.01 ≤ p < 0.05: %d (%.2f%%)\\n", bin3+0, 100*(bin3+0)/total
                printf "0.05 ≤ p < 0.1: %d (%.2f%%)\\n", bin4+0, 100*(bin4+0)/total
                printf "p ≥ 0.1: %d (%.2f%%)\\n", bin5+0, 100*(bin5+0)/total
            }
        }' >> ${prefix}.hwe_summary.txt
    fi
    
    echo "" >> ${prefix}.hwe_summary.txt
    
    # Worst violations
    if [ \$violations -gt 0 ]; then
        echo "TOP 10 HWE VIOLATIONS (lowest p-values):" >> ${prefix}.hwe_summary.txt
        echo "Chromosome\\tPosition\\tVariant_ID\\tREF\\tALT\\tHWE_P_Value" >> ${prefix}.hwe_summary.txt
        tail -n +2 ${prefix}.hwe_violations.txt | sort -k6,6g | head -10 >> ${prefix}.hwe_summary.txt
    fi
    
    echo "" >> ${prefix}.hwe_summary.txt
    echo "Detailed statistics: ${prefix}.hwe_stats.txt" >> ${prefix}.hwe_summary.txt
    echo "All violations: ${prefix}.hwe_violations.txt" >> ${prefix}.hwe_summary.txt
    
    # Add recommendation
    if [ \$violations -gt 0 ]; then
        echo "" >> ${prefix}.hwe_summary.txt
        echo "RECOMMENDATION:" >> ${prefix}.hwe_summary.txt
        echo "Consider filtering variants with severe HWE violations (p < ${pval_threshold})" >> ${prefix}.hwe_summary.txt
        echo "Review population stratification and sample quality" >> ${prefix}.hwe_summary.txt
    fi
    
    # Cleanup
    rm -f temp_hwe.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}