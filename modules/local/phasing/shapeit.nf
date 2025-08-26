process SHAPEIT_PHASE {
    tag "$meta.id"
    label 'process_high'
    
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path genetic_map
    
    output:
    tuple val(meta), path("${prefix}.phased.vcf.gz"), path("${prefix}.phased.vcf.gz.tbi"), emit: phased_vcf
    path "${prefix}.shapeit.log", emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def map_arg = genetic_map ? "--map $genetic_map" : ""
    """
    # Run SHAPEIT4 phasing
    shapeit4 \\
        --input $vcf \\
        --output ${prefix}.phased.vcf.gz \\
        --region ${meta.chr} \\
        $map_arg \\
        --thread $task.cpus \\
        $args \\
        2>&1 | tee ${prefix}.shapeit.log
    
    # Index the output
    bcftools index -t ${prefix}.phased.vcf.gz
    
    # Extract phasing statistics
    echo "\\nPhasing Statistics:" >> ${prefix}.shapeit.log
    bcftools query -f '[%SAMPLE\\t%GT\\n]' ${prefix}.phased.vcf.gz | \\
        awk -F'\\t' '{
            if(\$2 ~ /\\|/) phased++
            else unphased++
            total++
        } END {
            print "Total genotypes: "total
            print "Phased: "phased" ("100*phased/total"%)"
            print "Unphased: "unphased" ("100*unphased/total"%)"
        }' >> ${prefix}.shapeit.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit4: \$(shapeit4 --version 2>&1 | head -n1 | sed 's/.*v//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}