process MINIMAC4_IMPUTE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(vcf), val(ref_name), path(m3vcf), path(ref_vcf)
    val window
    val min_ratio

    output:
    tuple val(meta), path("*.dose.vcf.gz"), path("*.dose.vcf.gz.tbi"), emit: imputed_vcf
    tuple val(meta), path("*.info.gz"), emit: info
    tuple val(meta), path("*.erate"), emit: erate
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region = meta.region ? "--region ${meta.region}" : ''
    """
    minimac4 \\
        --refHaps ${m3vcf} \\
        --haps ${vcf} \\
        --prefix ${prefix} \\
        --format GT,DS,GP,HDS \\
        --noPhoneHome \\
        --log \\
        --window ${window} \\
        --min-ratio ${min_ratio} \\
        --cpus ${task.cpus} \\
        ${region} \\
        ${args}

    # Index the output VCF
    bcftools index -t ${prefix}.dose.vcf.gz

    # Compress info file
    gzip ${prefix}.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version 2>&1 | grep -E '^Minimac4' | sed 's/Minimac4 //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dose.vcf.gz
    touch ${prefix}.dose.vcf.gz.tbi
    touch ${prefix}.info.gz
    touch ${prefix}.erate
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: 4.1.6
        bcftools: 1.20
    END_VERSIONS
    """
}