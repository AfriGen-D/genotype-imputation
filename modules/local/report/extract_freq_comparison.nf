process EXTRACT_FREQ_COMPARISON {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), val(ref_name), path(imputed_vcf), path(imputed_index), path(ref_vcf)
    
    output:
    tuple val(meta), val(ref_name), path("*.freq_comparison.tsv"), emit: frequencies
    path "versions.yml"                                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    # Extract allele frequencies from imputed VCF
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/MAF\\t%INFO/R2\\n' \\
        ${imputed_vcf} | \\
    awk 'BEGIN {OFS="\\t"} 
         {
            # Handle missing values
            if (\$6 == ".") \$6 = "0";
            if (\$7 == ".") \$7 = "0";
            if (\$8 == ".") \$8 = "0";
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8
         }' > imputed_freq.tmp
    
    # Extract allele frequencies from reference panel VCF
    if [[ -f "${ref_vcf}" ]] && [[ "${ref_vcf}" != "NO_FILE" ]]; then
        bcftools query \\
            -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\n' \\
            ${ref_vcf} | \\
        awk 'BEGIN {OFS="\\t"} 
             {
                # Handle missing values
                if (\$6 == ".") \$6 = "0";
                # Calculate MAF from AF
                af = \$6
                maf = (af > 0.5) ? 1 - af : af
                print \$1, \$2, \$3, \$4, \$5, af, maf
             }' > ref_freq.tmp
        
        # Merge the two frequency files
        echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tIMP_AF\\tIMP_MAF\\tIMP_R2\\tREF_AF\\tREF_MAF\\tAF_DIFF\\tMAF_DIFF" > ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
        
        # Join on CHROM, POS, REF, ALT
        join -t \$'\\t' -1 1 -2 1 \\
            <(awk '{print \$1":"\$2":"\$4":"\$5"\\t"\$0}' imputed_freq.tmp | sort -k1,1) \\
            <(awk '{print \$1":"\$2":"\$4":"\$5"\\t"\$6"\\t"\$7}' ref_freq.tmp | sort -k1,1) | \\
        awk -F'\\t' '{
            # Parse the joined key
            split(\$1, key, ":");
            chrom = key[1];
            pos = key[2];
            ref = key[3];
            alt = key[4];
            
            # Imputed values
            id = \$4;
            imp_af = \$7;
            imp_maf = \$8;
            imp_r2 = \$9;
            
            # Reference values
            ref_af = \$10;
            ref_maf = \$11;
            
            # Calculate differences
            af_diff = imp_af - ref_af;
            maf_diff = imp_maf - ref_maf;
            
            print chrom, pos, id, ref, alt, imp_af, imp_maf, imp_r2, ref_af, ref_maf, af_diff, maf_diff
        }' OFS='\\t' >> ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
    else
        # No reference panel available, just output imputed frequencies
        echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tIMP_AF\\tIMP_MAF\\tIMP_R2\\tREF_AF\\tREF_MAF\\tAF_DIFF\\tMAF_DIFF" > ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
        awk '{print \$0"\\tNA\\tNA\\tNA\\tNA"}' imputed_freq.tmp >> ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
    fi
    
    # Get variant count for logging
    n_variants=\$(tail -n +2 ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv | wc -l)
    echo "Extracted and compared frequencies for \${n_variants} variants"
    
    # Clean up temp files
    rm -f imputed_freq.tmp ref_freq.tmp
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    touch ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.11
    END_VERSIONS
    """
}