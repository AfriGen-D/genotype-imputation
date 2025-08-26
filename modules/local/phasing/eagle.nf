process EAGLE_PHASE {
    tag "$meta.id:$meta.chr:$meta.region"
    label 'process_high'

    input:
    tuple val(meta), path(vcf), val(ref_name), path(ref_vcf)
    path genetic_map

    output:
    tuple val(meta), path("${prefix}.phased.vcf.gz"), path("${prefix}.phased.vcf.gz.tbi"), emit: phased_vcf
    tuple val(meta), path("${prefix}.eagle.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.chr}_${meta.region}"
    def region_arg = meta.region ? "--region ${meta.chr}:${meta.region}" : ""
    def map_arg = genetic_map ? "--geneticMapFile ${genetic_map}" : ""
    """
    # Run Eagle phasing
    eagle \\
        --vcfTarget ${vcf} \\
        --vcfRef ${ref_vcf} \\
        ${map_arg} \\
        --vcfOutFormat z \\
        --outPrefix ${prefix}.phased \\
        --numThreads ${task.cpus} \\
        --allowRefAltSwap \\
        --pbwtIters 2 \\
        ${region_arg} \\
        ${args} \\
        2>&1 | tee ${prefix}.eagle.log
    
    # Rename output (Eagle adds .vcf.gz automatically)
    if [ -f "${prefix}.phased.vcf.gz" ]; then
        echo "Output file exists"
    else
        mv ${prefix}.phased.bcf ${prefix}.phased.vcf.gz 2>/dev/null || true
    fi
    
    # Index output
    bcftools index -t ${prefix}.phased.vcf.gz
    
    # Extract phasing statistics
    echo "\\nPhasing Statistics:" >> ${prefix}.eagle.log
    bcftools query -f '[%SAMPLE\\t%GT\\n]' ${prefix}.phased.vcf.gz | \\
        awk -F'\\t' '{
            if(\$2 ~ /\\|/) phased++
            else unphased++
            total++
        } END {
            print "Total genotypes: "total
            print "Phased: "phased" ("100*phased/total"%)"
            print "Unphased: "unphased" ("100*unphased/total"%)"
        }' >> ${prefix}.eagle.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle: \$(eagle --version 2>&1 | grep -E '^Eagle' | sed 's/Eagle v//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}