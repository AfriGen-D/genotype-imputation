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
    
    # Calculate allele frequencies (handle missing INFO fields gracefully)
    # First check what INFO fields are available
    available_fields=\$(bcftools view -h ${vcf} | grep "^##INFO" | sed 's/.*ID=\\([^,]*\\).*/\\1/' | tr '\\n' ' ')
    echo "Available INFO fields: \$available_fields" >&2
    
    # Build query string based on available fields
    query_string='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT'
    if echo "\$available_fields" | grep -q "AF"; then
        query_string="\${query_string}\\t%INFO/AF"
    else
        query_string="\${query_string}\\t."
    fi
    if echo "\$available_fields" | grep -q "AC"; then
        query_string="\${query_string}\\t%INFO/AC"
    else
        query_string="\${query_string}\\t."
    fi
    if echo "\$available_fields" | grep -q "AN"; then
        query_string="\${query_string}\\t%INFO/AN"
    else
        query_string="\${query_string}\\t."
    fi
    query_string="\${query_string}\\n"
    
    bcftools query -f "\$query_string" ${vcf} > ${prefix}.frq
    
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
    
    # Apply basic QC filters based on available fields
    # Build filter expression based on what's available
    filter_expr=""
    
    # Check if QUAL field has actual values (not all missing)
    # Sample a few variants to check if QUAL is present
    has_qual=\$(bcftools query -f '%QUAL\\n' ${vcf} | head -100 | grep -v "^\\.\$" | wc -l)
    if [ "\$has_qual" -gt 0 ]; then
        echo "Found QUAL scores in VCF - will filter by QUAL>=20" >&2
        filter_expr="QUAL>=20"
    else
        echo "No QUAL scores found (all missing) - skipping QUAL filter" >&2
    fi
    
    # Check if DP exists in INFO
    if echo "\$available_fields" | grep -q "DP"; then
        if [ -n "\$filter_expr" ]; then
            filter_expr="\${filter_expr} && INFO/DP>=10"
        else
            filter_expr="INFO/DP>=10"
        fi
    fi
    
    # If no filters can be applied, just copy the input
    if [ -z "\$filter_expr" ]; then
        echo "No standard QC fields (QUAL, DP) found - applying site missingness filter only" >&2
        # Calculate and filter by site missingness if needed
        bcftools +fill-tags ${vcf} -Oz -o ${prefix}.qc.vcf.gz -- -t F_MISSING
    else
        echo "Applying filter: \$filter_expr" >&2
        bcftools filter \\
            --include "\$filter_expr" \\
            --output-type z \\
            --output ${prefix}.qc.vcf.gz \\
            --threads ${task.cpus} \\
            ${vcf}
    fi
    
    # Index the output
    bcftools index -t ${prefix}.qc.vcf.gz
    
    # Report QC summary
    echo "QC Summary for ${prefix}:" 
    input_count=\$(bcftools view -H ${vcf} | wc -l)
    output_count=\$(bcftools view -H ${prefix}.qc.vcf.gz | wc -l)
    echo "Total variants before QC: \$input_count"
    echo "Total variants after QC: \$output_count"
    echo "Variants filtered out: \$((input_count - output_count))"
    if [ -n "\$filter_expr" ]; then
        echo "Filters applied: \$filter_expr"
    else
        echo "Filters applied: Site missingness only (F_MISSING)"
    fi
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