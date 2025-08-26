process BEAGLE_PHASE {
    tag "$meta.id"
    label 'process_high'
    
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path genetic_map
    
    output:
    tuple val(meta), path("${prefix}.phased.vcf.gz"), path("${prefix}.phased.vcf.gz.tbi"), emit: phased_vcf
    path "${prefix}.beagle.log", emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def map_arg = genetic_map ? "map=$genetic_map" : ""
    """
    # Run Beagle phasing
    java -Xmx${task.memory.toGiga()}g -jar \${BEAGLE_JAR:-/usr/local/bin/beagle.jar} \\
        gt=$vcf \\
        out=${prefix}.phased \\
        chrom=${meta.chr} \\
        $map_arg \\
        nthreads=$task.cpus \\
        $args \\
        2>&1 | tee ${prefix}.beagle.log
    
    # Beagle outputs .vcf.gz by default
    if [ -f ${prefix}.phased.vcf.gz ]; then
        mv ${prefix}.phased.vcf.gz ${prefix}.phased.vcf.gz.tmp
        bcftools view ${prefix}.phased.vcf.gz.tmp -Oz -o ${prefix}.phased.vcf.gz
        rm ${prefix}.phased.vcf.gz.tmp
    fi
    
    # Index the output
    bcftools index -t ${prefix}.phased.vcf.gz
    
    # Extract phasing statistics
    echo "\\nPhasing Statistics:" >> ${prefix}.beagle.log
    bcftools query -f '[%SAMPLE\\t%GT\\n]' ${prefix}.phased.vcf.gz | \\
        awk -F'\\t' '{
            if(\$2 ~ /\\|/) phased++
            else unphased++
            total++
        } END {
            print "Total genotypes: "total
            print "Phased: "phased" ("100*phased/total"%)"
            print "Unphased: "unphased" ("100*unphased/total"%)"
        }' >> ${prefix}.beagle.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: 5.4
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}