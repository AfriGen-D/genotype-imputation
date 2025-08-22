process VERIFY_GENETIC_MAP {
    tag "genetic_map_verification"
    label 'process_single'
    
    container 'quay.io/biocontainers/bcftools:1.11--h7c999a4_0'
    
    input:
    path genetic_map
    path input_vcf
    path reference_genome  // Optional - for build detection
    
    output:
    path "genetic_map_verification.txt", emit: report
    path "versions.yml"                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Detect genome build from filenames and paths
    echo "==================================" > genetic_map_verification.txt
    echo "GENOME BUILD DETECTION" >> genetic_map_verification.txt
    echo "==================================" >> genetic_map_verification.txt
    
    # Detect genetic map build
    genetic_map_build="unknown"
    if echo "${genetic_map}" | grep -q "hg19\\|b37\\|GRCh37"; then
        genetic_map_build="hg19/b37"
    elif echo "${genetic_map}" | grep -q "hg38\\|b38\\|GRCh38"; then
        genetic_map_build="hg38/b38"
    fi
    echo "Genetic map appears to be: \$genetic_map_build" >> genetic_map_verification.txt
    
    # Detect reference genome build if provided
    ref_build="not_provided"
    if [[ "${reference_genome}" != "NO_FILE" ]] && [[ -n "${reference_genome}" ]]; then
        if echo "${reference_genome}" | grep -q "hg19\\|b37\\|GRCh37"; then
            ref_build="hg19/b37"
        elif echo "${reference_genome}" | grep -q "hg38\\|b38\\|GRCh38\\|assembly38"; then
            ref_build="hg38/b38"
        else
            ref_build="unknown"
        fi
        echo "Reference genome appears to be: \$ref_build" >> genetic_map_verification.txt
    fi
    
    # Detect input data build from filename
    input_build="unknown"
    if echo "${input_vcf}" | grep -q "b37\\|hg19\\|GRCh37"; then
        input_build="hg19/b37"
    elif echo "${input_vcf}" | grep -q "b38\\|hg38\\|GRCh38"; then
        input_build="hg38/b38"
    fi
    echo "Input VCF appears to be: \$input_build (based on filename)" >> genetic_map_verification.txt
    
    # Check for build consistency
    echo "" >> genetic_map_verification.txt
    if [[ "\$genetic_map_build" != "unknown" ]] && [[ "\$ref_build" != "unknown" ]] && [[ "\$ref_build" != "not_provided" ]]; then
        if [[ "\$genetic_map_build" != "\$ref_build" ]]; then
            echo "WARNING: Genetic map build (\$genetic_map_build) does not match reference genome build (\$ref_build)" >> genetic_map_verification.txt
        fi
    fi
    
    if [[ "\$genetic_map_build" != "unknown" ]] && [[ "\$input_build" != "unknown" ]]; then
        if [[ "\$genetic_map_build" != "\$input_build" ]]; then
            echo "WARNING: Genetic map build (\$genetic_map_build) may not match input data build (\$input_build)" >> genetic_map_verification.txt
        fi
    fi
    
    echo "" >> genetic_map_verification.txt
    echo "==================================" >> genetic_map_verification.txt
    echo "GENETIC MAP COVERAGE CHECK" >> genetic_map_verification.txt
    echo "==================================" >> genetic_map_verification.txt
    
    # Get unique chromosomes and their position ranges from VCF
    bcftools query -f '%CHROM\\t%POS\\n' ${input_vcf} | \\
    awk '{
        if (!seen[\$1]) {
            chr[\$1] = \$1
            min_pos[\$1] = \$2
            max_pos[\$1] = \$2
            seen[\$1] = 1
        } else {
            if (\$2 < min_pos[\$1]) min_pos[\$1] = \$2
            if (\$2 > max_pos[\$1]) max_pos[\$1] = \$2
        }
    } END {
        for (c in chr) {
            print c, min_pos[c], max_pos[c]
        }
    }' | sort -V > vcf_ranges.txt
    
    # Check genetic map coverage for each chromosome
    echo "Checking genetic map coverage..." >> genetic_map_verification.txt
    echo "=================================" >> genetic_map_verification.txt
    
    # Decompress genetic map if needed
    if [[ ${genetic_map} == *.gz ]]; then
        zcat ${genetic_map} > genetic_map.txt
    else
        cp ${genetic_map} genetic_map.txt
    fi
    
    # Check each chromosome's coverage
    while read -r chr min_pos max_pos; do
        echo "" >> genetic_map_verification.txt
        echo "Chromosome: \$chr" >> genetic_map_verification.txt
        echo "  VCF range: \$min_pos - \$max_pos" >> genetic_map_verification.txt
        
        # Try both with and without 'chr' prefix
        chr_noprefix=\${chr#chr}
        chr_withprefix="chr\$chr_noprefix"
        
        # Check genetic map coverage
        map_min=\$(awk -v c1="\$chr" -v c2="\$chr_noprefix" -v c3="\$chr_withprefix" \\
            '(\$1==c1 || \$1==c2 || \$1==c3) {print \$2}' genetic_map.txt | head -1)
        map_max=\$(awk -v c1="\$chr" -v c2="\$chr_noprefix" -v c3="\$chr_withprefix" \\
            '(\$1==c1 || \$1==c2 || \$1==c3) {print \$2}' genetic_map.txt | tail -1)
        
        if [[ -z "\$map_min" ]] || [[ -z "\$map_max" ]]; then
            echo "  WARNING: No genetic map data found for chromosome \$chr" >> genetic_map_verification.txt
            echo "  Status: MISSING" >> genetic_map_verification.txt
        else
            echo "  Genetic map range: \$map_min - \$map_max" >> genetic_map_verification.txt
            
            # Check if VCF range is covered
            if [[ \$min_pos -ge \$map_min ]] && [[ \$max_pos -le \$map_max ]]; then
                echo "  Status: OK (fully covered)" >> genetic_map_verification.txt
            elif [[ \$max_pos -lt \$map_min ]] || [[ \$min_pos -gt \$map_max ]]; then
                echo "  Status: ERROR (no overlap)" >> genetic_map_verification.txt
                echo "  ERROR: VCF positions fall completely outside genetic map coverage!" >> genetic_map_verification.txt
            else
                # Calculate coverage percentage
                overlap_start=\$(( min_pos > map_min ? min_pos : map_min ))
                overlap_end=\$(( max_pos < map_max ? max_pos : map_max ))
                
                if [[ \$overlap_end -ge \$overlap_start ]]; then
                    vcf_span=\$(( max_pos - min_pos + 1 ))
                    overlap_span=\$(( overlap_end - overlap_start + 1 ))
                    coverage_pct=\$(( overlap_span * 100 / vcf_span ))
                    echo "  Status: WARNING (partial coverage: \${coverage_pct}%)" >> genetic_map_verification.txt
                    
                    if [[ \$min_pos -lt \$map_min ]]; then
                        echo "  WARNING: VCF starts at \$min_pos but genetic map starts at \$map_min" >> genetic_map_verification.txt
                    fi
                    if [[ \$max_pos -gt \$map_max ]]; then
                        echo "  WARNING: VCF ends at \$max_pos but genetic map ends at \$map_max" >> genetic_map_verification.txt
                        echo "  This will cause Eagle phasing errors for positions beyond \$map_max" >> genetic_map_verification.txt
                    fi
                else
                    echo "  Status: ERROR (no valid overlap)" >> genetic_map_verification.txt
                fi
            fi
        fi
    done < vcf_ranges.txt
    
    # Add summary
    echo "" >> genetic_map_verification.txt
    echo "=================================" >> genetic_map_verification.txt
    echo "SUMMARY" >> genetic_map_verification.txt
    echo "=================================" >> genetic_map_verification.txt
    
    # Check for any errors or warnings
    if grep -q "ERROR" genetic_map_verification.txt; then
        echo "CRITICAL: Genetic map coverage errors detected!" >> genetic_map_verification.txt
        echo "Eagle phasing will fail for regions without genetic map coverage." >> genetic_map_verification.txt
        echo "" >> genetic_map_verification.txt
        echo "SOLUTION: Use a genetic map that matches your data's genome build:" >> genetic_map_verification.txt
        echo "  - For hg19/b37 data: use genetic_map_hg19_withX.txt.gz" >> genetic_map_verification.txt
        echo "  - For hg38/b38 data: use genetic_map_hg38_withX.txt.gz" >> genetic_map_verification.txt
        
        # Exit with error to stop pipeline
        echo "" >> genetic_map_verification.txt
        echo "Pipeline stopped due to genetic map incompatibility." >> genetic_map_verification.txt
        cat genetic_map_verification.txt >&2
        exit 1
    elif grep -q "WARNING" genetic_map_verification.txt; then
        echo "Warnings detected - some regions may have limited phasing accuracy." >> genetic_map_verification.txt
    else
        echo "All chromosomes have full genetic map coverage - OK to proceed." >> genetic_map_verification.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    echo "Genetic map verification (stub mode)" > genetic_map_verification.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.11
    END_VERSIONS
    """
}