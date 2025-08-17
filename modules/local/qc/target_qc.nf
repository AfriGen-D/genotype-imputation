process TARGET_QC {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.qc.vcf.gz"), path("*.qc.vcf.gz.tbi"), emit: vcf
    path "*.stats"                                                , emit: stats
    path "*.frq"                                                  , emit: freq
    path "*.missing"                                              , emit: missing
    path "versions.yml"                                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Calculate comprehensive statistics
    bcftools stats \\
        --threads ${task.cpus} \\
        ${vcf} > ${prefix}.stats
    
    # Calculate allele frequencies
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/AC\\t%INFO/AN\\n' \\
        ${vcf} > ${prefix}.frq
    
    # Calculate missing data per variant
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t[%GT\\t]\\n' \\
        ${vcf} | \\
        awk '{
            total = 0
            missing = 0
            for(i=4; i<=NF; i++) {
                total++
                if(\$i == "./.") missing++
            }
            if(total > 0) {
                miss_rate = missing/total
                print \$1"\\t"\$2"\\t"\$3"\\t"missing"\\t"total"\\t"miss_rate
            }
        }' > ${prefix}.missing
    
    # Apply basic QC filters and create output VCF
    bcftools filter \\
        --include 'QUAL>=20 && INFO/DP>=10' \\
        --output-type z \\
        --output ${prefix}.qc.vcf.gz \\
        --threads ${task.cpus} \\
        ${vcf}
    
    # Index the output
    bcftools index -t ${prefix}.qc.vcf.gz
    
    # Report QC summary
    echo "QC Summary for ${prefix}:" 
    echo "Total variants before QC: \$(bcftools view -H ${vcf} | wc -l)"
    echo "Total variants after QC: \$(bcftools view -H ${prefix}.qc.vcf.gz | wc -l)"
    echo "Average missing rate: \$(awk '{sum+=\$6; count++} END {if(count>0) print sum/count; else print 0}' ${prefix}.missing)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qc.vcf.gz
    touch ${prefix}.qc.vcf.gz.tbi
    touch ${prefix}.stats
    touch ${prefix}.frq
    touch ${prefix}.missing
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}