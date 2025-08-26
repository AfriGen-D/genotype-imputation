process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_medium'
    
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path annotation_db
    
    output:
    tuple val(meta), path("${prefix}.annotated.vcf.gz"), path("${prefix}.annotated.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.annotation_stats.txt"), emit: stats
    tuple val(meta), path("${prefix}.annotation_summary.txt"), emit: report
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Get initial variant counts
    n_input=\$(bcftools view -H ${vcf} | wc -l)
    echo "Input variants: \$n_input" > ${prefix}.annotation_stats.txt
    
    # Add basic variant annotations using bcftools +fill-tags
    bcftools +fill-tags ${vcf} -Oz -o temp_tagged.vcf.gz -- -t all
    bcftools index -t temp_tagged.vcf.gz
    
    # Annotate with external database if provided
    if [ -f "${annotation_db}" ]; then
        echo "Annotating with external database: ${annotation_db}" >> ${prefix}.annotation_stats.txt
        
        # Check if annotation database is VCF or other format
        if [[ "${annotation_db}" == *.vcf.gz ]] || [[ "${annotation_db}" == *.vcf ]]; then
            # VCF-based annotation
            bcftools annotate \\
                --annotations ${annotation_db} \\
                --columns ID,INFO \\
                --output-type z \\
                --output ${prefix}.annotated.vcf.gz \\
                --threads ${task.cpus} \\
                ${args} \\
                temp_tagged.vcf.gz
        else
            # Assume tab-delimited format
            bcftools annotate \\
                --annotations ${annotation_db} \\
                --header-lines <(echo '##INFO=<ID=ANNOTATION,Number=1,Type=String,Description="External annotation">') \\
                --columns CHROM,POS,REF,ALT,ANNOTATION \\
                --output-type z \\
                --output ${prefix}.annotated.vcf.gz \\
                --threads ${task.cpus} \\
                ${args} \\
                temp_tagged.vcf.gz
        fi
    else
        echo "No external annotation database provided - using basic annotations only" >> ${prefix}.annotation_stats.txt
        cp temp_tagged.vcf.gz ${prefix}.annotated.vcf.gz
    fi
    
    # Index output
    bcftools index -t ${prefix}.annotated.vcf.gz
    
    # Get annotation statistics
    n_output=\$(bcftools view -H ${prefix}.annotated.vcf.gz | wc -l)
    echo "Output variants: \$n_output" >> ${prefix}.annotation_stats.txt
    
    # Count variants with annotations
    n_with_id=\$(bcftools query -f '%ID\\n' ${prefix}.annotated.vcf.gz | grep -v '^[.]\$' | wc -l)
    n_with_qual=\$(bcftools query -f '%QUAL\\n' ${prefix}.annotated.vcf.gz | grep -v '^[.]\$' | wc -l)
    
    echo "Variants with ID: \$n_with_id" >> ${prefix}.annotation_stats.txt
    echo "Variants with QUAL: \$n_with_qual" >> ${prefix}.annotation_stats.txt
    
    # Check for specific INFO tags
    for tag in AC AN AF MAF HWE; do
        if bcftools view -h ${prefix}.annotated.vcf.gz | grep -q "ID=\$tag"; then
            n_with_tag=\$(bcftools query -i "\$tag!='.' && \$tag!=''" -f '%CHROM\\n' ${prefix}.annotated.vcf.gz 2>/dev/null | wc -l)
            echo "Variants with \$tag: \$n_with_tag" >> ${prefix}.annotation_stats.txt
        fi
    done
    
    # Generate summary report
    echo "Variant Annotation Summary for ${meta.id}" > ${prefix}.annotation_summary.txt
    echo "Generated on: \$(date)" >> ${prefix}.annotation_summary.txt
    echo "" >> ${prefix}.annotation_summary.txt
    
    echo "INPUT/OUTPUT SUMMARY:" >> ${prefix}.annotation_summary.txt
    cat ${prefix}.annotation_stats.txt >> ${prefix}.annotation_summary.txt
    echo "" >> ${prefix}.annotation_summary.txt
    
    # Annotation coverage rates
    if [ \$n_output -gt 0 ]; then
        echo "ANNOTATION COVERAGE:" >> ${prefix}.annotation_summary.txt
        id_rate=\$(echo "scale=2; 100 * \$n_with_id / \$n_output" | bc)
        qual_rate=\$(echo "scale=2; 100 * \$n_with_qual / \$n_output" | bc)
        echo "Variants with ID: \$id_rate%" >> ${prefix}.annotation_summary.txt
        echo "Variants with QUAL: \$qual_rate%" >> ${prefix}.annotation_summary.txt
        
        for tag in AC AN AF MAF HWE; do
            if bcftools view -h ${prefix}.annotated.vcf.gz | grep -q "ID=\$tag"; then
                n_with_tag=\$(bcftools query -i "\$tag!='.' && \$tag!=''" -f '%CHROM\\n' ${prefix}.annotated.vcf.gz 2>/dev/null | wc -l)
                tag_rate=\$(echo "scale=2; 100 * \$n_with_tag / \$n_output" | bc)
                echo "Variants with \$tag: \$tag_rate%" >> ${prefix}.annotation_summary.txt
            fi
        done
    fi
    
    echo "" >> ${prefix}.annotation_summary.txt
    
    # INFO fields summary
    echo "AVAILABLE INFO FIELDS:" >> ${prefix}.annotation_summary.txt
    bcftools view -h ${prefix}.annotated.vcf.gz | grep '^##INFO' | cut -d',' -f1 | sed 's/.*ID=//' >> ${prefix}.annotation_summary.txt
    
    echo "" >> ${prefix}.annotation_summary.txt
    
    # Sample annotation stats
    n_samples=\$(bcftools query -l ${prefix}.annotated.vcf.gz | wc -l)
    echo "SAMPLE INFORMATION:" >> ${prefix}.annotation_summary.txt
    echo "Number of samples: \$n_samples" >> ${prefix}.annotation_summary.txt
    
    # FORMAT fields
    echo "" >> ${prefix}.annotation_summary.txt
    echo "AVAILABLE FORMAT FIELDS:" >> ${prefix}.annotation_summary.txt
    bcftools view -h ${prefix}.annotated.vcf.gz | grep '^##FORMAT' | cut -d',' -f1 | sed 's/.*ID=//' >> ${prefix}.annotation_summary.txt
    
    echo "" >> ${prefix}.annotation_summary.txt
    echo "Output file: ${prefix}.annotated.vcf.gz" >> ${prefix}.annotation_summary.txt
    echo "Detailed statistics: ${prefix}.annotation_stats.txt" >> ${prefix}.annotation_summary.txt
    
    # Cleanup
    rm -f temp_tagged.vcf.gz temp_tagged.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}