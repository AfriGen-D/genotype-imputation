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
    path "*.phasing_strategy.txt"                                         , emit: strategy, optional: true
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
    
    # Try Eagle phasing with reference panel BCF
    set +e  # Don't exit on error immediately
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
    
    EAGLE_EXIT_CODE=\$?
    set -e
    
    # Check if Eagle failed due to genetic map issues
    if [ \$EAGLE_EXIT_CODE -ne 0 ]; then
        if grep -q "Genetic distance range.*0 cM" ${prefix}_${chunk_id}.phasing.log || \\
           grep -q "genetic distance ranges must be positive" ${prefix}_${chunk_id}.phasing.log; then
            echo "WARNING: Eagle failed due to zero genetic distance in region ${chunk_id}"
            echo "This chunk (${chrm}:${chunk_start}-${chunk_end}) appears to be in a low-recombination region"
            echo "Attempting alternative strategies to preserve this region..."
            
            # Strategy 1: Try phasing without region boundaries (use whole chromosome context)
            echo "Strategy 1: Phasing without region boundaries to use broader genetic context..."
            set +e
            eagle \\
                --vcfTarget=${vcf} \\
                --geneticMapFile=${genetic_map} \\
                ${has_ref ? "--vcfRef=${reference_panel}" : ""} \\
                --vcfOutFormat=z \\
                --noImpMissing \\
                --numThreads=${task.cpus} \\
                --pbwtIters=${params.eagle_pbwt_iters ?: 2} \\
                --chrom=${chrm} \\
                --outPrefix=${prefix}_${chunk_id}.phased \\
                $args \\
                2>&1 | tee ${prefix}_${chunk_id}.phasing_retry1.log
            
            RETRY1_EXIT_CODE=\$?
            set -e
            
            if [ \$RETRY1_EXIT_CODE -eq 0 ]; then
                echo "SUCCESS: Phased using whole chromosome context"
                echo "STRATEGY: whole_chromosome" > ${prefix}_${chunk_id}.phasing_strategy.txt
                echo "REASON: zero_genetic_distance" >> ${prefix}_${chunk_id}.phasing_strategy.txt
            else
                # Strategy 2: Try with extended flanking regions (double the buffer)
                echo "Strategy 2: Expanding flanking regions for more genetic context..."
                extended_buffer=\$((${params.buffer_size ?: 1000000} * 2))
                extended_start=\$((${chunk_start} - \$extended_buffer))
                extended_end=\$((${chunk_end} + \$extended_buffer))
                
                # Ensure start is not negative
                if [ \$extended_start -lt 1 ]; then
                    extended_start=1
                fi
                
                set +e
                eagle \\
                    --vcfTarget=${vcf} \\
                    --geneticMapFile=${genetic_map} \\
                    ${has_ref ? "--vcfRef=${reference_panel}" : ""} \\
                    --vcfOutFormat=z \\
                    --noImpMissing \\
                    --numThreads=${task.cpus} \\
                    --pbwtIters=${params.eagle_pbwt_iters ?: 2} \\
                    --chrom=${chrm} \\
                    --bpStart=\$extended_start \\
                    --bpEnd=\$extended_end \\
                    --bpFlanking=\$extended_buffer \\
                    --outPrefix=${prefix}_${chunk_id}.phased \\
                    $args \\
                    2>&1 | tee ${prefix}_${chunk_id}.phasing_retry2.log
                
                RETRY2_EXIT_CODE=\$?
                set -e
                
                if [ \$RETRY2_EXIT_CODE -eq 0 ]; then
                    echo "SUCCESS: Phased using extended flanking regions"
                    echo "STRATEGY: extended_flanking" > ${prefix}_${chunk_id}.phasing_strategy.txt
                    echo "REASON: zero_genetic_distance" >> ${prefix}_${chunk_id}.phasing_strategy.txt
                    echo "BUFFER: \$extended_buffer" >> ${prefix}_${chunk_id}.phasing_strategy.txt
                else
                    # Strategy 3: Try without genetic map (uses LD patterns only)
                    echo "Strategy 3: Phasing without genetic map (LD-based only)..."
                    set +e
                    eagle \\
                        --vcfTarget=${vcf} \\
                        ${has_ref ? "--vcfRef=${reference_panel}" : ""} \\
                        --vcfOutFormat=z \\
                        --noImpMissing \\
                        --numThreads=${task.cpus} \\
                        --pbwtIters=${params.eagle_pbwt_iters ?: 2} \\
                        --chrom=${chrm} \\
                        --bpStart=${chunk_start} \\
                        --bpEnd=${chunk_end} \\
                        --bpFlanking=${params.buffer_size ?: 1000000} \\
                        --outPrefix=${prefix}_${chunk_id}.phased \\
                        $args \\
                        2>&1 | tee ${prefix}_${chunk_id}.phasing_retry3.log
                    
                    RETRY3_EXIT_CODE=\$?
                    set -e
                    
                    if [ \$RETRY3_EXIT_CODE -eq 0 ]; then
                        echo "SUCCESS: Phased using LD patterns without genetic map"
                        echo "WARNING: Phasing quality may be reduced in low-recombination region" >> ${prefix}_${chunk_id}.phasing.log
                        echo "STRATEGY: ld_only" > ${prefix}_${chunk_id}.phasing_strategy.txt
                        echo "REASON: zero_genetic_distance" >> ${prefix}_${chunk_id}.phasing_strategy.txt
                        echo "WARNING: no_genetic_map" >> ${prefix}_${chunk_id}.phasing_strategy.txt
                    else
                        echo "ERROR: All phasing strategies failed for this region"
                        echo "This region may require specialized handling or alternative phasing tools"
                        # Exit with special code to mark for further investigation
                        exit 199
                    fi
                fi
            fi
        else
            # Eagle failed for other reasons
            echo "ERROR: Eagle failed with exit code \$EAGLE_EXIT_CODE"
            exit \$EAGLE_EXIT_CODE
        fi
    fi
    
    # Create strategy file for normal phasing if successful
    if [ \$EAGLE_EXIT_CODE -eq 0 ]; then
        echo "STRATEGY: standard" > ${prefix}_${chunk_id}.phasing_strategy.txt
        echo "REASON: normal_phasing" >> ${prefix}_${chunk_id}.phasing_strategy.txt
    fi
    
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