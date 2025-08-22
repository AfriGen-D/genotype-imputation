process CHECK_BUILD_CONSISTENCY {
    tag "build_check"
    label 'process_single'
    
    container 'quay.io/biocontainers/bcftools:1.11--h7c999a4_0'
    
    input:
    path genetic_map
    path reference_genome
    
    output:
    path "build_consistency.txt", emit: report
    path "versions.yml"         , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    echo "==================================" > build_consistency.txt
    echo "GENOME BUILD CONSISTENCY CHECK" >> build_consistency.txt
    echo "==================================" >> build_consistency.txt
    echo "" >> build_consistency.txt
    
    # Detect genetic map build
    genetic_map_build="unknown"
    if echo "${genetic_map}" | grep -q "hg19\\|b37\\|GRCh37"; then
        genetic_map_build="hg19/b37"
    elif echo "${genetic_map}" | grep -q "hg38\\|b38\\|GRCh38"; then
        genetic_map_build="hg38/b38"
    fi
    echo "Genetic map: ${genetic_map}" >> build_consistency.txt
    echo "  Detected build: \$genetic_map_build" >> build_consistency.txt
    echo "" >> build_consistency.txt
    
    # Detect reference genome build
    ref_build="unknown"
    if echo "${reference_genome}" | grep -q "hg19\\|b37\\|GRCh37"; then
        ref_build="hg19/b37"
    elif echo "${reference_genome}" | grep -q "hg38\\|b38\\|GRCh38\\|assembly38"; then
        ref_build="hg38/b38"
    fi
    echo "Reference genome: ${reference_genome}" >> build_consistency.txt
    echo "  Detected build: \$ref_build" >> build_consistency.txt
    echo "" >> build_consistency.txt
    
    # Check consistency
    echo "==================================" >> build_consistency.txt
    echo "CONSISTENCY CHECK RESULT" >> build_consistency.txt
    echo "==================================" >> build_consistency.txt
    
    if [[ "\$genetic_map_build" == "unknown" ]]; then
        echo "WARNING: Could not detect genetic map build from filename" >> build_consistency.txt
        echo "  Please verify manually that genetic map matches your data" >> build_consistency.txt
    fi
    
    if [[ "\$ref_build" == "unknown" ]]; then
        echo "WARNING: Could not detect reference genome build from filename" >> build_consistency.txt
        echo "  Please verify manually that reference genome matches your data" >> build_consistency.txt
    fi
    
    if [[ "\$genetic_map_build" != "unknown" ]] && [[ "\$ref_build" != "unknown" ]]; then
        if [[ "\$genetic_map_build" == "\$ref_build" ]]; then
            echo "✓ PASS: Genetic map and reference genome builds match (\$genetic_map_build)" >> build_consistency.txt
        else
            echo "✗ FAIL: Build mismatch detected!" >> build_consistency.txt
            echo "" >> build_consistency.txt
            echo "  Genetic map build: \$genetic_map_build" >> build_consistency.txt
            echo "  Reference genome build: \$ref_build" >> build_consistency.txt
            echo "" >> build_consistency.txt
            echo "This mismatch will cause errors during phasing and imputation." >> build_consistency.txt
            echo "" >> build_consistency.txt
            echo "SOLUTION:" >> build_consistency.txt
            echo "  For hg38/b38 data, use:" >> build_consistency.txt
            echo "    eagle_genetic_map = '/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s/map/genetic_map_hg38_withX.txt.gz'" >> build_consistency.txt
            echo "  For hg19/b37 data, use:" >> build_consistency.txt
            echo "    eagle_genetic_map = '/cbio/dbs/refpanels/eagle/tables/genetic_map_hg19_withX.txt.gz'" >> build_consistency.txt
            
            cat build_consistency.txt >&2
            exit 1
        fi
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    echo "Build consistency check (stub mode)" > build_consistency.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.11
    END_VERSIONS
    """
}