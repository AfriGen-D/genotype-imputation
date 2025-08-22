process GET_REF_CHROMOSOMES {
    tag "reference_panel"
    label 'process_single'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    val ref_panel_info  // Array with [name, msav_template, bcf_template]
    val expected_build  // Expected genome build
    
    output:
    path "ref_chromosomes.txt", emit: chromosomes
    path "ref_build.txt"      , emit: ref_build
    path "versions.yml"       , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def ref_name = ref_panel_info[0]
    def bcf_template = ref_panel_info[2]
    """
    # Get list of available reference panel chromosomes
    echo "Scanning reference panel: ${ref_name}" >&2
    
    # Check which chromosome files exist
    touch ref_chromosomes.txt
    detected_build=""
    
    # Check common chromosome names with 'chr' prefix (b38 style)
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \\
               chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \\
               chr20 chr21 chr22 chrX chrY chrM chrMT; do
        ref_file="${bcf_template.replace('%s', "\${chr}")}"
        if [ -f "\${ref_file}" ]; then
            echo "\${chr}" >> ref_chromosomes.txt
            echo "  Found: \${chr}" >&2
            detected_build="b38"
        fi
    done
    
    # Check without 'chr' prefix (b37 style)
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT; do
        ref_file="${bcf_template.replace('%s', "\${chr}")}"
        if [ -f "\${ref_file}" ]; then
            echo "\${chr}" >> ref_chromosomes.txt
            echo "  Found: \${chr}" >&2
            if [ -z "\${detected_build}" ]; then
                detected_build="b37"
            fi
        fi
    done
    
    # Write detected build
    echo "\${detected_build}" > ref_build.txt
    echo "Reference panel build detected: \${detected_build}" >&2
    
    # Check if reference panel build matches expected build
    if [ "${expected_build}" != "auto" ] && [ "\${detected_build}" != "${expected_build}" ]; then
        echo "================================================================" >&2
        echo "ERROR: Reference panel build mismatch!" >&2
        echo "  Reference panel appears to be: \${detected_build}" >&2
        echo "  Expected build from config: ${expected_build}" >&2
        echo "" >&2
        echo "Please ensure your reference panel matches the expected build:" >&2
        echo "  - b37: chromosomes as '1', '2', etc." >&2
        echo "  - b38: chromosomes as 'chr1', 'chr2', etc." >&2
        echo "================================================================" >&2
        exit 1
    fi
    
    # Report summary
    echo "Total chromosomes found: \$(wc -l < ref_chromosomes.txt)" >&2
    
    if [ ! -s ref_chromosomes.txt ]; then
        echo "ERROR: No reference panel files found matching pattern: ${bcf_template}" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
    
    stub:
    """
    echo "chr1" > ref_chromosomes.txt
    echo "chr2" >> ref_chromosomes.txt
    echo "chr21" >> ref_chromosomes.txt
    echo "chr22" >> ref_chromosomes.txt
    echo "chrX" >> ref_chromosomes.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.0
    END_VERSIONS
    """
}