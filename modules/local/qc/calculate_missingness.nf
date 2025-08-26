process CALCULATE_MISSINGNESS {
    tag "$meta.id"
    label 'process_single'
    
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("${prefix}.missingness.txt"), emit: report
    tuple val(meta), path("${prefix}.sample_missingness.txt"), emit: sample_stats
    tuple val(meta), path("${prefix}.variant_missingness.txt"), emit: variant_stats
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Calculate per-sample missingness
    bcftools query -f '[%SAMPLE\\t%GT\\n]' ${vcf} | \\
        awk '
        BEGIN { OFS="\\t" }
        {
            samples[\$1]++
            if(\$2 == "./." || \$2 == ".|.") missing[\$1]++
        }
        END {
            print "Sample\\tTotal_Sites\\tMissing_Sites\\tMissingness_Rate" > "${prefix}.sample_missingness.txt"
            for(sample in samples) {
                miss = missing[sample] + 0
                rate = samples[sample] > 0 ? miss/samples[sample] : 0
                print sample "\\t" samples[sample] "\\t" miss "\\t" rate >> "${prefix}.sample_missingness.txt"
            }
        }'
    
    # Calculate per-variant missingness
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t[%GT\\t]\\n' ${vcf} | \\
        awk '
        BEGIN { OFS="\\t" }
        {
            variant = \$1":"\$2
            variant_id = \$3 != "." ? \$3 : variant
            total = 0
            missing = 0
            for(i=4; i<=NF; i++) {
                if(\$i != "") {
                    total++
                    if(\$i == "./." || \$i == ".|.") missing++
                }
            }
            if(total > 0) {
                rate = missing/total
                print \$1 "\\t" \$2 "\\t" variant_id "\\t" total "\\t" missing "\\t" rate
            }
        }' > temp_variant_miss.txt
    
    # Add header and sort by missingness rate
    echo -e "Chromosome\\tPosition\\tVariant_ID\\tTotal_Samples\\tMissing_Samples\\tMissingness_Rate" > ${prefix}.variant_missingness.txt
    sort -k6,6nr temp_variant_miss.txt >> ${prefix}.variant_missingness.txt
    
    # Generate summary report
    echo "Missingness Analysis for ${meta.id}" > ${prefix}.missingness.txt
    echo "Generated on: \$(date)" >> ${prefix}.missingness.txt
    echo "" >> ${prefix}.missingness.txt
    
    # Sample missingness summary
    echo "SAMPLE MISSINGNESS SUMMARY:" >> ${prefix}.missingness.txt
    total_samples=\$(tail -n +2 ${prefix}.sample_missingness.txt | wc -l)
    echo "Total samples: \$total_samples" >> ${prefix}.missingness.txt
    
    if [ \$total_samples -gt 0 ]; then
        avg_sample_miss=\$(tail -n +2 ${prefix}.sample_missingness.txt | awk '{sum+=\$4} END {print sum/NR}')
        max_sample_miss=\$(tail -n +2 ${prefix}.sample_missingness.txt | awk 'NR==1{max=\$4} {if(\$4>max) max=\$4} END {print max}')
        high_miss_samples=\$(tail -n +2 ${prefix}.sample_missingness.txt | awk '\$4 > 0.1 {count++} END {print count+0}')
        
        echo "Average sample missingness: \$avg_sample_miss" >> ${prefix}.missingness.txt
        echo "Maximum sample missingness: \$max_sample_miss" >> ${prefix}.missingness.txt
        echo "Samples with >10% missingness: \$high_miss_samples" >> ${prefix}.missingness.txt
    fi
    
    echo "" >> ${prefix}.missingness.txt
    
    # Variant missingness summary  
    echo "VARIANT MISSINGNESS SUMMARY:" >> ${prefix}.missingness.txt
    total_variants=\$(tail -n +2 ${prefix}.variant_missingness.txt | wc -l)
    echo "Total variants: \$total_variants" >> ${prefix}.missingness.txt
    
    if [ \$total_variants -gt 0 ]; then
        avg_variant_miss=\$(tail -n +2 ${prefix}.variant_missingness.txt | awk '{sum+=\$6} END {print sum/NR}')
        max_variant_miss=\$(tail -n +2 ${prefix}.variant_missingness.txt | awk 'NR==1{max=\$6} {if(\$6>max) max=\$6} END {print max}')
        high_miss_variants=\$(tail -n +2 ${prefix}.variant_missingness.txt | awk '\$6 > 0.1 {count++} END {print count+0}')
        
        echo "Average variant missingness: \$avg_variant_miss" >> ${prefix}.missingness.txt
        echo "Maximum variant missingness: \$max_variant_miss" >> ${prefix}.missingness.txt
        echo "Variants with >10% missingness: \$high_miss_variants" >> ${prefix}.missingness.txt
    fi
    
    echo "" >> ${prefix}.missingness.txt
    echo "Detailed per-sample statistics: ${prefix}.sample_missingness.txt" >> ${prefix}.missingness.txt
    echo "Detailed per-variant statistics: ${prefix}.variant_missingness.txt" >> ${prefix}.missingness.txt
    
    # Cleanup
    rm -f temp_variant_miss.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}