process IMPUTE5_IMPUTE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(vcf), val(ref_name), path(ref_vcf), path(ref_map)
    val buffer_size
    val ne
    val min_ratio

    output:
    tuple val(meta), path("*.imputed.vcf.gz"), path("*.imputed.vcf.gz.tbi"), emit: imputed_vcf
    tuple val(meta), path("*.info"), emit: info
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region = meta.region ? "--region ${meta.region}" : ''
    def map_arg = ref_map ? "--map ${ref_map}" : ''
    """
    impute5 \\
        --h ${ref_vcf} \\
        --g ${vcf} \\
        --r ${meta.chr}:${meta.start}-${meta.end} \\
        --o ${prefix}.imputed.vcf.gz \\
        --l ${prefix}.log \\
        --buffer-size ${buffer_size} \\
        --ne ${ne} \\
        --min-ratio ${min_ratio} \\
        --threads ${task.cpus} \\
        ${map_arg} \\
        ${region} \\
        ${args}

    # Index the output VCF
    bcftools index -t ${prefix}.imputed.vcf.gz

    # Extract info scores
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/INFO\\n' \\
        ${prefix}.imputed.vcf.gz > ${prefix}.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        impute5: \$(impute5 --version 2>&1 | grep -E '^IMPUTE5' | sed 's/IMPUTE5 //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.imputed.vcf.gz
    touch ${prefix}.imputed.vcf.gz.tbi
    touch ${prefix}.info
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        impute5: 1.2.0
        bcftools: 1.20
    END_VERSIONS
    """
}