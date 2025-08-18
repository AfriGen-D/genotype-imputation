process EXTRACT_FREQ_REFERENCE {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), val(ref_name), path(ref_vcf)
    
    output:
    tuple val(meta), val(ref_name), path("*.ref_freq.tsv"), emit: frequencies
    path "versions.yml"                                    , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    # Check if reference VCF exists and is valid
    if [[ -f "${ref_vcf}" ]] && [[ "${ref_vcf}" != "NO_FILE" ]]; then
        # Extract allele frequencies from reference panel VCF
        bcftools query \\
            -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\n' \\
            ${ref_vcf} | \\
        awk 'BEGIN {
                OFS="\\t"
                print "CHROM", "POS", "ID", "REF", "ALT", "REF_AF", "REF_MAF"
            } 
            {
                # Handle missing values
                if (\$6 == ".") \$6 = "0";
                
                # Calculate MAF from AF
                af = \$6
                maf = (af > 0.5) ? 1 - af : af
                
                # Output with REF_ prefix for clarity
                print \$1, \$2, \$3, \$4, \$5, af, maf
            }' > ${prefix}_${ref_name}_${chunk_id}.ref_freq.tsv
        
        # Get variant count for logging
        n_variants=\$(tail -n +2 ${prefix}_${ref_name}_${chunk_id}.ref_freq.tsv | wc -l)
        echo "Extracted frequencies for \${n_variants} variants from reference panel"
    else
        # Create empty file with header if no reference panel available
        echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tREF_AF\\tREF_MAF" > ${prefix}_${ref_name}_${chunk_id}.ref_freq.tsv
        echo "No reference panel VCF available, created empty frequency file"
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tREF_AF\\tREF_MAF" > ${prefix}_${ref_name}_${chunk_id}.ref_freq.tsv
    echo -e "chr21\\t10000001\\trs123\\tA\\tG\\t0.04\\t0.04" >> ${prefix}_${ref_name}_${chunk_id}.ref_freq.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}