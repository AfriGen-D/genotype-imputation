process EXTRACT_FREQ_SIMPLE {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), val(ref_name), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), val(ref_name), path("*.freq.tsv"), emit: frequencies
    path "versions.yml"                                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    # Extract allele frequencies from imputed VCF
    # Note: This extracts from imputed data only - reference panel comparison to be added
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/MAF\\t%INFO/R2\\n' \\
        ${vcf} | \\
    awk 'BEGIN {OFS="\\t"; print "CHROM","POS","ID","REF","ALT","AF","MAF","R2"} 
         {
            # Handle missing values
            if (\$6 == ".") \$6 = "0";
            if (\$7 == ".") \$7 = "0";
            if (\$8 == ".") \$8 = "0";
            print \$0
         }' > ${prefix}_${ref_name}_${chunk_id}.freq.tsv
    
    # Get variant count for logging
    n_variants=\$(tail -n +2 ${prefix}_${ref_name}_${chunk_id}.freq.tsv | wc -l)
    echo "Extracted frequencies for \${n_variants} variants from imputed data"
    
    # TODO: Add reference panel frequency extraction when ref VCFs are available
    echo "Note: Reference panel frequency comparison not yet implemented"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    touch ${prefix}_${ref_name}_${chunk_id}.freq.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.11
    END_VERSIONS
    """
}