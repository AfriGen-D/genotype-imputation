process QC_SITE_MISSINGNESS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.nomiss.vcf.gz"), path("*.nomiss.vcf.gz.tbi"), emit: vcf
    path "*.missingness.txt"                                              , emit: missingness
    path "versions.yml"                                                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_missing = params.site_missingness ?: 0.05
    """
    # First check what INFO fields are available
    available_fields=\$(bcftools view -h ${vcf} | grep "^##INFO" | sed 's/.*ID=\\([^,]*\\).*/\\1/' | tr '\\n' ' ')
    echo "Available INFO fields: \$available_fields" >&2
    
    # Build query string based on available fields
    query_string='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT'
    if echo "\$available_fields" | grep -q "AN"; then
        query_string="\${query_string}\\t%INFO/AN"
    else
        query_string="\${query_string}\\t."
    fi
    if echo "\$available_fields" | grep -q "AC"; then
        query_string="\${query_string}\\t%INFO/AC"
    else
        query_string="\${query_string}\\t."
    fi
    query_string="\${query_string}\\n"
    
    # Extract site information (for documentation purposes)
    bcftools query -f "\$query_string" ${vcf} > ${prefix}.info.txt
    
    # First ensure F_MISSING tag is calculated if not present
    if ! bcftools view -h ${vcf} | grep -q "##INFO=<ID=F_MISSING"; then
        echo "Adding F_MISSING tag to VCF" >&2
        bcftools +fill-tags ${vcf} -Oz -o ${prefix}.tagged.vcf.gz -- -t F_MISSING
        bcftools index -t ${prefix}.tagged.vcf.gz
        input_vcf="${prefix}.tagged.vcf.gz"
    else
        input_vcf="${vcf}"
    fi
    
    # Filter sites with high missingness
    bcftools filter \\
        --include "F_MISSING<${max_missing}" \\
        --output-type z \\
        --output ${prefix}.nomiss.vcf.gz \\
        \$input_vcf
    
    # Index the output
    bcftools index -t ${prefix}.nomiss.vcf.gz
    
    # Report missingness statistics
    echo "Site missingness statistics:" > ${prefix}.missingness.txt
    echo "Threshold: ${max_missing}" >> ${prefix}.missingness.txt
    echo "Sites before filtering:" >> ${prefix}.missingness.txt
    bcftools view -H \$input_vcf | wc -l >> ${prefix}.missingness.txt
    echo "Sites after filtering:" >> ${prefix}.missingness.txt
    bcftools view -H ${prefix}.nomiss.vcf.gz | wc -l >> ${prefix}.missingness.txt
    
    # Calculate average missingness if possible
    if bcftools view -h ${prefix}.nomiss.vcf.gz | grep -q "##INFO=<ID=F_MISSING"; then
        echo "Average site missingness after filtering:" >> ${prefix}.missingness.txt
        bcftools query -f '%INFO/F_MISSING\\n' ${prefix}.nomiss.vcf.gz | \\
            awk '{sum+=\$1; count++} END {if(count>0) print sum/count; else print "N/A"}' >> ${prefix}.missingness.txt
    fi
    
    # Clean up temporary files
    if [ -f "${prefix}.tagged.vcf.gz" ]; then
        rm -f ${prefix}.tagged.vcf.gz ${prefix}.tagged.vcf.gz.tbi
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nomiss.vcf.gz
    touch ${prefix}.nomiss.vcf.gz.tbi
    touch ${prefix}.missingness.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}