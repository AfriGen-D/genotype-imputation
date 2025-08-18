process EXTRACT_FREQ_IMPUTED {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), val(ref_name), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), val(ref_name), path("*.imputed_freq.tsv"), emit: frequencies
    path "versions.yml"                                        , emit: versions
    
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
        ${vcf} | \\
    awk 'BEGIN {
            OFS="\\t"
            print "CHROM", "POS", "ID", "REF", "ALT", "IMP_AF", "IMP_MAF", "IMP_R2"
        } 
        {
            # Handle missing values
            if (\$6 == ".") \$6 = "0";
            if (\$7 == ".") \$7 = "0";
            if (\$8 == ".") \$8 = "0";
            
            # Output with IMP_ prefix for clarity
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8
        }' > ${prefix}_${ref_name}_${chunk_id}.imputed_freq.tsv
    
    # Get variant count for logging
    n_variants=\$(tail -n +2 ${prefix}_${ref_name}_${chunk_id}.imputed_freq.tsv | wc -l)
    echo "Extracted frequencies for \${n_variants} variants from imputed data"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tIMP_AF\\tIMP_MAF\\tIMP_R2" > ${prefix}_${ref_name}_${chunk_id}.imputed_freq.tsv
    echo -e "chr21\\t10000001\\trs123\\tA\\tG\\t0.05\\t0.05\\t0.95" >> ${prefix}_${ref_name}_${chunk_id}.imputed_freq.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}