process VALIDATE_IMPUTATION_CHUNK {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.1.0'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(ref_name), path(ref_msav), path(ref_vcf)
    
    output:
    tuple val(meta), path("*.validation.txt"), emit: report
    tuple val(meta), path("*.validation.status"), emit: status
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer_size = params.buffer_size ?: 500000
    def min_ratio = params.minRatio ?: 0.01
    
    // Extract chromosome and coordinates
    if (!meta.contig) error "ERROR: meta.contig is required for VALIDATE_IMPUTATION_CHUNK"
    def chrom = meta.contig
    def start = meta.start ?: error("ERROR: meta.start is required")
    def end = meta.end ?: error("ERROR: meta.end is required")
    
    """
    # Validate the chunk for imputation suitability
    python3 /users/mamana/genotype-imputation/bin/validate_imputation_chunk.py \\
        --vcf ${vcf} \\
        --ref-vcf ${ref_vcf} \\
        --chrom ${chrom} \\
        --start ${start} \\
        --end ${end} \\
        --buffer ${buffer_size} \\
        --min-ratio ${min_ratio} \\
        --output ${prefix}.validation.txt \\
        --status-file ${prefix}.validation.status \\
        ${args}
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "PASS" > ${prefix}.validation.status
    touch ${prefix}.validation.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        bcftools: 1.20
    END_VERSIONS
    """
}