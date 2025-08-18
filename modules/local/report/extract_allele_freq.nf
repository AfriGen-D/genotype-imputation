process EXTRACT_ALLELE_FREQ {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), val(ref_name), path(vcf), path(vcf_index), val(source)
    
    output:
    tuple val(meta), val(ref_name), val(source), path("*.freq.tsv"), emit: frequencies
    path "versions.yml"                                             , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    def chrm = meta.contig ?: ''
    """
    # Extract allele frequencies from VCF
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/MAF\\n' \\
        ${vcf} | \\
    awk 'BEGIN {OFS="\\t"; print "CHROM","POS","ID","REF","ALT","AF","MAF","SOURCE"} 
         {
            # Handle missing values
            if (\$6 == ".") \$6 = "0";
            if (\$7 == ".") \$7 = "0";
            print \$0,"${source}"
         }' > ${prefix}_${ref_name}_${chunk_id}_${source}.freq.tsv
    
    # Get variant count for logging
    n_variants=\$(tail -n +2 ${prefix}_${ref_name}_${chunk_id}_${source}.freq.tsv | wc -l)
    echo "Extracted frequencies for \${n_variants} variants from ${source} data"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    touch ${prefix}_${ref_name}_${chunk_id}_${source}.freq.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.11
    END_VERSIONS
    """
}