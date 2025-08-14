process BEAGLE5_IMPUTE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(vcf), val(ref_name), path(ref_vcf)
    path genetic_map
    val ne
    val window

    output:
    tuple val(meta), path("*.imputed.vcf.gz"), path("*.imputed.vcf.gz.tbi"), emit: imputed_vcf
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def map_arg = genetic_map ? "map=${genetic_map}" : ''
    def region = meta.region ? "chrom=${meta.chr}:${meta.start}-${meta.end}" : ''
    """
    java -Xmx${task.memory.toGiga()}g -jar \$BEAGLE_JAR \\
        gt=${vcf} \\
        ref=${ref_vcf} \\
        out=${prefix} \\
        ne=${ne} \\
        window=${window} \\
        nthreads=${task.cpus} \\
        ${map_arg} \\
        ${region} \\
        ${args} \\
        2>&1 | tee ${prefix}.log

    # Rename output (Beagle adds .vcf.gz automatically)
    mv ${prefix}.vcf.gz ${prefix}.imputed.vcf.gz

    # Index the output VCF
    bcftools index -t ${prefix}.imputed.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(java -jar \$BEAGLE_JAR 2>&1 | grep -E '^Beagle' | sed 's/Beagle //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.imputed.vcf.gz
    touch ${prefix}.imputed.vcf.gz.tbi
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: 5.4
        bcftools: 1.20
    END_VERSIONS
    """
}