process IMPUTE_MINIMAC4 {
    tag "$meta.id"
    label 'process_high'
    
    container 'mamana/imputation:minimac4-4.1.6'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(ref_name), path(ref_msav), path(ref_vcf)
    
    output:
    tuple val(meta), val(ref_name), path("*.dose.vcf.gz"), path("*.dose.vcf.gz.tbi"), emit: imputed
    tuple val(meta), path("*.sites.vcf.gz")                                          , emit: info
    path "versions.yml"                                                               , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    def chrm = meta.contig ?: ''
    def chunk_start = meta.start ?: ''
    def chunk_end = meta.end ?: ''
    
    // Build region string if we have chunk coordinates
    def region_arg = ''
    if (chrm && chunk_start && chunk_end) {
        region_arg = "--region ${chrm}:${chunk_start}-${chunk_end}"
    }
    
    """
    # Run Minimac4 imputation
    minimac4 \\
        --refHaps $ref_msav \\
        --haps $vcf \\
        --prefix ${prefix}_${ref_name}_${chunk_id} \\
        --format GT,DS,GP \\
        --noPhoneHome \\
        --minRatio ${params.minRatio} \\
        --cpus $task.cpus \\
        ${region_arg} \\
        $args
    
    # Index the output VCF
    bcftools index -t ${prefix}_${ref_name}_${chunk_id}.dose.vcf.gz
    
    # Create version file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version 2>&1 | head -1 | sed 's/^.*Minimac4 - //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    touch ${prefix}_${ref_name}_${chunk_id}.dose.vcf.gz
    touch ${prefix}_${ref_name}_${chunk_id}.dose.vcf.gz.tbi
    touch ${prefix}_${ref_name}_${chunk_id}.sites.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: 4.1.6
        bcftools: 1.20
    END_VERSIONS
    """
}