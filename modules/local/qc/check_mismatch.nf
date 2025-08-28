process CHECK_MISMATCH {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), path(chunk_vcf), path(reference_panel), path(reference_panel_index)
    
    output:
    tuple val(meta), path("*.mismatch.txt"), path("*.mismatch.status"), emit: mismatch
    path "versions.yml"                                                , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!meta.contig) error "ERROR: meta.contig is required for CHECK_MISMATCH but was not provided"
    def chrm = meta.contig
    def chunk_start = meta.start ?: ''
    def chunk_end = meta.end ?: ''
    def max_mismatch_rate = params.max_mismatch_rate ?: 0.2  // Maximum 20% mismatch allowed (relaxed for diverse populations)
    def min_matching_variants = params.min_matching_variants ?: 20  // Minimum matching variants (relaxed from 50)
    """
    # Extract variants from chunk
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${chunk_vcf} > ${prefix}.chunk.alleles
    
    # Extract variants from reference panel in chunk region
    if [ -n "${chunk_start}" ] && [ -n "${chunk_end}" ]; then
        bcftools query -r ${chrm}:${chunk_start}-${chunk_end} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${reference_panel} > ${prefix}.ref.alleles
    else
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${reference_panel} > ${prefix}.ref.alleles
    fi
    
    # Check allele concordance using awk
    awk '
    BEGIN {
        total_chunk = 0
        total_ref = 0
        matching_pos = 0
        matching_alleles = 0
        mismatching_alleles = 0
        chunk_only = 0
        ref_only = 0
    }
    # Read reference alleles
    NR==FNR {
        ref[\$1"\\t"\$2] = \$3"\\t"\$4
        total_ref++
        next
    }
    # Process chunk alleles
    {
        total_chunk++
        key = \$1"\\t"\$2
        if (key in ref) {
            matching_pos++
            if (ref[key] == \$3"\\t"\$4) {
                matching_alleles++
            } else {
                mismatching_alleles++
                print "Mismatch at", key":", \$3, \$4, "vs", ref[key] > "${prefix}.mismatches.detail"
            }
        } else {
            chunk_only++
        }
    }
    END {
        for (key in ref) {
            if (!(key in processed)) {
                ref_only++
            }
        }
        
        print "Chunk variants:", total_chunk
        print "Reference variants in region:", total_ref
        print "Positions in both:", matching_pos
        print "Matching alleles:", matching_alleles
        print "Mismatching alleles:", mismatching_alleles
        print "Chunk-only variants:", chunk_only
        print "Reference-only variants:", ref_only
        
        if (matching_pos > 0) {
            mismatch_rate = mismatching_alleles / matching_pos
            print "Mismatch rate:", mismatch_rate
            
            # Write status based on mismatch rate with relaxed criteria
            if (mismatch_rate <= ${max_mismatch_rate} && matching_alleles >= ${min_matching_variants}) {
                print "PASS" > "${prefix}.mismatch.status"
                print "Status: PASS - Acceptable mismatch rate and sufficient matching variants"
            } else if (mismatch_rate <= 0.3 && matching_alleles >= 10) {
                # Relaxed criteria for challenging regions
                print "WARN" > "${prefix}.mismatch.status"
                print "Status: WARN - Marginal match quality, proceeding with caution"
                print "Warning: Mismatch rate", mismatch_rate, "with", matching_alleles, "matching variants"
            } else {
                print "FAIL" > "${prefix}.mismatch.status"
                if (mismatch_rate > ${max_mismatch_rate}) {
                    print "Status: FAIL - High mismatch rate:", mismatch_rate
                } else {
                    print "Status: FAIL - Too few matching variants:", matching_alleles
                }
            }
        } else {
            print "FAIL" > "${prefix}.mismatch.status"
            print "Status: FAIL - No overlapping positions"
        }
    }
    ' ${prefix}.ref.alleles ${prefix}.chunk.alleles > ${prefix}.mismatch.txt
    
    # Clean up temporary files
    rm -f ${prefix}.chunk.alleles ${prefix}.ref.alleles
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mismatch.txt
    echo "PASS" > ${prefix}.mismatch.status
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}