process CHECK_GENOME_BUILD {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), path(vcf)
    val expected_build  // 'b37' or 'b38' from params
    
    output:
    tuple val(meta), path(vcf), emit: vcf
    path "*.build_check.txt"   , emit: build_status
    path "versions.yml"        , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Get the first chromosome from the VCF to check naming convention
    # Use head -n1 and ignore SIGPIPE (exit code 141) from bcftools
    first_chrom=\$(bcftools query -f '%CHROM\\n' ${vcf} 2>/dev/null | head -n1 || true)
    
    # Make sure we got a chromosome
    if [ -z "\$first_chrom" ]; then
        echo "ERROR: Could not read chromosome from VCF file ${vcf}" >&2
        exit 1
    fi
    
    echo "Dataset: ${meta.id}" > ${prefix}.build_check.txt
    echo "First chromosome found: \$first_chrom" >> ${prefix}.build_check.txt
    echo "Expected genome build: ${expected_build}" >> ${prefix}.build_check.txt
    
    # Determine the build based on chromosome naming
    if [[ "\$first_chrom" =~ ^chr ]]; then
        detected_build="b38"
        echo "Detected build: b38 (chromosomes have 'chr' prefix)" >> ${prefix}.build_check.txt
    else
        detected_build="b37"
        echo "Detected build: b37 (chromosomes without 'chr' prefix)" >> ${prefix}.build_check.txt
    fi
    
    # Check if detected build matches expected build
    if [ "${expected_build}" == "auto" ]; then
        echo "Build check: AUTO mode - accepting detected build \$detected_build" >> ${prefix}.build_check.txt
        echo "PASS" >> ${prefix}.build_check.txt
    elif [ "\$detected_build" == "${expected_build}" ]; then
        echo "Build check: PASS - detected build matches expected build" >> ${prefix}.build_check.txt
        echo "PASS" >> ${prefix}.build_check.txt
    else
        echo "================================================================" >&2
        echo "WARNING: Genome build mismatch detected!" >&2
        echo "Dataset: ${meta.id}" >&2
        echo "  Detected build: \$detected_build (based on chromosome naming)" >&2
        echo "  Expected build: ${expected_build} (from configuration)" >&2
        echo "" >&2
        echo "Genome build conventions:" >&2
        echo "  - b37/hg19: chromosomes as '1', '2', '3', etc." >&2
        echo "  - b38/hg38: chromosomes as 'chr1', 'chr2', 'chr3', etc." >&2
        echo "" >&2
        echo "This dataset will be EXCLUDED from the analysis." >&2
        echo "================================================================" >&2
        
        echo "Build check: FAIL - build mismatch" >> ${prefix}.build_check.txt
        echo "FAIL" >> ${prefix}.build_check.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "PASS" > ${prefix}.build_check.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}