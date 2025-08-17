process EAGLE_PHASING {
    tag "$meta.id"
    label 'process_high'
    
    container 'mamana/eagle-vcf-processing:eagle-2.4.1'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path genetic_map
    path reference_panel  
    path reference_panel_index    
    output:
    tuple val(meta), path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi"), emit: phased
    path "*.log"                                                          , emit: log
    path "versions.yml"                                                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    def chrm = meta.contig ?: ''
    def chunk_start = meta.start ?: ''
    def chunk_end = meta.end ?: ''
    
    // Check if we have reference panel and genetic map for reference-based phasing
    def has_ref = reference_panel && reference_panel.name != 'NO_FILE'
    def has_map = genetic_map && genetic_map.name != 'NO_FILE'
    
    // Genetic map is REQUIRED for Eagle phasing
    if (!has_map) {
        error "ERROR: Eagle genetic map file is required but not provided. Please set eagle_genetic_map parameter in your config file."
    }
    
    // Build Eagle command arguments
    def eagle_args = []
    if (chrm) {
        eagle_args << "--chrom=${chrm}"
    }
    if (chunk_start && chunk_end) {
        eagle_args << "--bpStart=${chunk_start}"
        eagle_args << "--bpEnd=${chunk_end}"
        eagle_args << "--bpFlanking=${params.buffer_size ?: 1000000}"
    }
    def eagle_args_str = eagle_args.join(' ')
    
    """
    # Index the input VCF if needed
    if [ ! -f ${vcf}.tbi ]; then
        bcftools index -t ${vcf}
    fi
    
    # Run Eagle phasing with reference panel BCF
    eagle \\
        --vcfTarget=${vcf} \\
        --geneticMapFile=${genetic_map} \\
        ${has_ref ? "--vcfRef=${reference_panel}" : ""} \\
        --vcfOutFormat=z \\
        --noImpMissing \\
        --numThreads=${task.cpus} \\
        --pbwtIters=${params.eagle_pbwt_iters ?: 2} \\
        ${eagle_args_str} \\
        --outPrefix=${prefix}_${chunk_id}.phased \\
        $args \\
        2>&1 | tee ${prefix}_${chunk_id}.phasing.log
    
    # Index the output VCF
    bcftools index -t ${prefix}_${chunk_id}.phased.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle: \$(eagle --version 2>&1 | head -1 | sed 's/^.*Eagle //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    touch ${prefix}_${chunk_id}.phased.vcf.gz
    touch ${prefix}_${chunk_id}.phased.vcf.gz.tbi
    touch ${prefix}_${chunk_id}.phasing.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle: 2.4.1
        bcftools: 1.20
    END_VERSIONS
    """
}