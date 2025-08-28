process IMPUTE_MINIMAC4_SAFE {
    tag "$meta.id"
    label 'process_high'
    
    container 'mamana/imputation:minimac4-4.1.6'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index), path(overlap_file, stageAs: "overlap.txt")
    tuple val(ref_name), path(ref_msav), path(ref_vcf)
    
    output:
    tuple val(meta), val(ref_name), path("*.dose.vcf.gz"), path("*.dose.vcf.gz.tbi"), optional: true, emit: imputed
    tuple val(meta), path("*.sites.vcf.gz"), optional: true, emit: info
    tuple val(meta), path("*.skipped.txt"), optional: true, emit: skipped
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    def min_ratio = params.minRatio ?: 0.01
    def min_variants = 50  // Minimum variants needed for imputation
    
    // Chromosome is required for imputation
    if (!meta.contig) error "ERROR: meta.contig is required for IMPUTE_MINIMAC4_SAFE but was not provided"
    def chrm = meta.contig
    def chunk_start = meta.start ?: ''
    def chunk_end = meta.end ?: ''
    
    """
    # First, check if this chunk has sufficient variants for imputation
    echo "Checking variant density for chunk ${chrm}:${chunk_start}-${chunk_end}..."
    
    # Count variants in the target VCF for this region
    VARIANT_COUNT=\$(bcftools view -H -r ${chrm}:${chunk_start}-${chunk_end} ${vcf} | wc -l)
    echo "Found \${VARIANT_COUNT} variants in target chunk"
    
    # Handle BCF index - create if missing
    if [ ! -f "${ref_vcf}.csi" ] && [ ! -f "\$(basename ${ref_vcf}).csi" ]; then
        echo "Creating BCF index..."
        bcftools index ${ref_vcf}
    fi
    
    # Count variants in reference panel for this region  
    REF_COUNT=\$(bcftools view -H -r ${chrm}:${chunk_start}-${chunk_end} ${ref_vcf} | wc -l)
    echo "Found \${REF_COUNT} variants in reference panel"
    
    # Calculate ratio
    if [ "\${REF_COUNT}" -gt 0 ]; then
        RATIO=\$(echo "scale=4; \${VARIANT_COUNT} / \${REF_COUNT}" | bc -l)
    else
        RATIO=0
    fi
    echo "Variant ratio: \${RATIO}"
    
    # Check if we should proceed with imputation
    MIN_RATIO="${min_ratio}"
    MIN_VARIANTS="${min_variants}"
    
    # Use awk for floating point comparison
    PROCEED=\$(awk -v ratio="\${RATIO}" -v min_ratio="\${MIN_RATIO}" -v var_count="\${VARIANT_COUNT}" -v min_var="\${MIN_VARIANTS}" '
        BEGIN {
            if (ratio >= min_ratio && var_count >= min_var) {
                print "YES"
            } else {
                print "NO"
            }
        }
    ')
    
    if [ "\${PROCEED}" = "YES" ]; then
        echo "Proceeding with imputation (ratio: \${RATIO}, variants: \${VARIANT_COUNT})"
        
        # Build region argument - use only the main chunk, not extended buffers
        # This prevents Minimac4 from trying to impute adjacent empty regions
        REGION_ARG="--region ${chrm}:${chunk_start}-${chunk_end}"
        
        # Run Minimac4 with error handling
        set +e  # Don't exit immediately on error
        
        minimac4 \\
            --refHaps ${ref_msav} \\
            --haps ${vcf} \\
            --prefix ${prefix}_${ref_name}_${chunk_id} \\
            --format GT,DS,GP \\
            --noPhoneHome \\
            --minRatio ${min_ratio} \\
            --cpus ${task.cpus} \\
            \${REGION_ARG} \\
            ${args} 2>&1 | tee minimac4.log
        
        MINIMAC_EXIT=\${PIPESTATUS[0]}
        
        if [ "\${MINIMAC_EXIT}" -eq 0 ]; then
            # Success - index the output
            echo "Imputation completed successfully"
            bcftools index -t ${prefix}_${ref_name}_${chunk_id}.dose.vcf.gz
        else
            # Check if it failed due to empty regions
            if grep -q "not enough target variants" minimac4.log; then
                echo "Imputation failed due to insufficient variants - creating skip marker"
                echo "SKIPPED: Not enough variants for imputation" > ${prefix}_${ref_name}_${chunk_id}.skipped.txt
                echo "Chunk: ${chrm}:${chunk_start}-${chunk_end}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
                echo "Target variants: \${VARIANT_COUNT}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
                echo "Reference variants: \${REF_COUNT}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
                echo "Ratio: \${RATIO}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
            else
                # Other error - let it fail normally
                echo "Imputation failed with unexpected error"
                exit \${MINIMAC_EXIT}
            fi
        fi
        
        set -e  # Re-enable exit on error
    else
        echo "Skipping imputation - insufficient variants"
        echo "  Ratio (\${RATIO}) < minimum (${min_ratio}) or"
        echo "  Variants (\${VARIANT_COUNT}) < minimum (${min_variants})"
        
        # Create skip marker file
        echo "SKIPPED: Insufficient variants for imputation" > ${prefix}_${ref_name}_${chunk_id}.skipped.txt
        echo "Chunk: ${chrm}:${chunk_start}-${chunk_end}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
        echo "Target variants: \${VARIANT_COUNT}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
        echo "Reference variants: \${REF_COUNT}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
        echo "Ratio: \${RATIO}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
        echo "Minimum ratio: ${min_ratio}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
        echo "Minimum variants: ${min_variants}" >> ${prefix}_${ref_name}_${chunk_id}.skipped.txt
    fi
    
    # Create version file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version 2>&1 | head -1 | sed 's/^.*Minimac4 - //' || echo "4.1.6")
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    touch ${prefix}_${ref_name}_${chunk_id}.dose.vcf.gz
    touch ${prefix}_${ref_name}_${chunk_id}.dose.vcf.gz.tbi
    touch ${prefix}_${ref_name}_${chunk_id}.sites.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: 4.1.6
        bcftools: 1.20
    END_VERSIONS
    """
}