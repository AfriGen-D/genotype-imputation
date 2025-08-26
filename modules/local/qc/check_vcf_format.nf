process CHECK_VCF_FORMAT {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.validated.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.validation.log"), emit: log
    tuple val(meta), path("${meta.id}.timing.log"), emit: timing
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def skip_stats = params.skip_vcf_stats ? 'true' : 'false'
    """
    export SKIP_STATS=${skip_stats}
    # Initialize timing log
    echo "=== CHECK_VCF_FORMAT Timing Report ===" > ${prefix}.timing.log
    echo "Process started at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${prefix}.timing.log
    START_TIME=\$(date +%s)
    
    # CRITICAL: Check for GT field (required for imputation)
    echo "Step 1: GT field check started at \$(date '+%H:%M:%S')" >> ${prefix}.timing.log
    STEP_START=\$(date +%s)
    echo "Validating VCF format..." > ${prefix}.validation.log
    # Check if GT field exists
    bcftools view -h ${vcf} | grep "##FORMAT=<ID=GT" > /dev/null 2>&1
    if [ \$? -ne 0 ]; then
        echo "ERROR: VCF missing GT format field" >> ${prefix}.validation.log
        exit 1
    fi
    echo "VCF has GT field - OK" >> ${prefix}.validation.log
    echo "  GT field check completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${prefix}.timing.log
    
    # CRITICAL: Create validated copy (just symlink if already compressed)
    echo "Step 2: Creating validated VCF started at \$(date '+%H:%M:%S')" >> ${prefix}.timing.log
    STEP_START=\$(date +%s)
    if [[ ${vcf} == *.gz ]]; then
        # Already compressed, just create symlink with absolute path
        ln -sf \$(readlink -f ${vcf}) ${prefix}.validated.vcf.gz
    else
        # Compress if not already
        bcftools view ${vcf} -Oz -o ${prefix}.validated.vcf.gz --threads ${task.cpus}
    fi
    echo "  Validated VCF created in \$((\$(date +%s) - STEP_START)) seconds" >> ${prefix}.timing.log
    
    # CRITICAL: Index the validated VCF (required for downstream)
    echo "Step 3: VCF indexing started at \$(date '+%H:%M:%S')" >> ${prefix}.timing.log
    STEP_START=\$(date +%s)
    # Check if index already exists for the original file
    ORIGINAL_VCF=\$(readlink -f ${vcf})
    if [[ -f "\${ORIGINAL_VCF}.tbi" ]] || [[ -f "\${ORIGINAL_VCF}.csi" ]]; then
        echo "  Using existing index from input file" >> ${prefix}.timing.log
        if [[ -f "\${ORIGINAL_VCF}.tbi" ]]; then
            ln -sf "\${ORIGINAL_VCF}.tbi" ${prefix}.validated.vcf.gz.tbi
        else
            ln -sf "\${ORIGINAL_VCF}.csi" ${prefix}.validated.vcf.gz.csi
        fi
    else
        echo "  Creating new index..." >> ${prefix}.timing.log
        bcftools index -t ${prefix}.validated.vcf.gz --threads ${task.cpus}
    fi
    echo "  Indexing completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${prefix}.timing.log
    
    # OPTIONAL: Quick stats (can be disabled if slow)
    if [ "\${SKIP_STATS:-false}" != "true" ]; then
        echo "Step 4: Quick stats started at \$(date '+%H:%M:%S')" >> ${prefix}.timing.log
        STEP_START=\$(date +%s)
        echo "\\nVCF Statistics:" >> ${prefix}.validation.log
        bcftools query -f 'x\\n' ${prefix}.validated.vcf.gz | wc -l | xargs echo "Total variants:" >> ${prefix}.validation.log
        echo "  Quick stats completed in \$((\$(date +%s) - STEP_START)) seconds" >> ${prefix}.timing.log
    fi
    
    # Final timing report
    END_TIME=\$(date +%s)
    echo "Process completed at: \$(date '+%Y-%m-%d %H:%M:%S')" >> ${prefix}.timing.log
    echo "Total execution time: \$((END_TIME - START_TIME)) seconds" >> ${prefix}.timing.log
    echo "=== End of Timing Report ===" >> ${prefix}.timing.log
    
    # Display timing summary
    cat ${prefix}.timing.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}