process CHECK_VCF_FORMAT {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.validated.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.validation.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Check VCF format and integrity
    bcftools view -h ${vcf} > header_check.txt
    
    # Validate VCF
    bcftools +fixploidy ${vcf} -- -d 2 2>&1 | tee ${prefix}.validation.log
    
    # Check for required fields
    if ! bcftools view -h ${vcf} | grep -q "##FORMAT=<ID=GT"; then
        echo "ERROR: VCF missing GT format field" >> ${prefix}.validation.log
        exit 1
    fi
    
    # Ensure proper compression and indexing
    bcftools view ${vcf} -Oz -o ${prefix}.validated.vcf.gz
    bcftools index -t ${prefix}.validated.vcf.gz
    
    # Get stats
    bcftools stats ${prefix}.validated.vcf.gz >> ${prefix}.validation.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}