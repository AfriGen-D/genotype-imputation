process MERGE_PHASED_CHUNKS {
    tag "$meta.id"
    label 'process_medium'
    
    
    input:
    tuple val(meta), path(vcf_files), path(tbi_files)
    
    output:
    tuple val(meta), path("${prefix}.merged.vcf.gz"), path("${prefix}.merged.vcf.gz.tbi"), emit: merged_vcf
    path "${prefix}.merge.log", emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Sort VCF files by chromosome and position (based on filename if they contain position info)
    ls *.vcf.gz | sort -V > vcf_list.txt
    
    # Count input files
    n_files=\$(cat vcf_list.txt | wc -l)
    echo "Merging \$n_files VCF chunks for ${meta.id}" > ${prefix}.merge.log
    echo "Input files:" >> ${prefix}.merge.log
    cat vcf_list.txt >> ${prefix}.merge.log
    
    if [ \$n_files -eq 1 ]; then
        # Single file, just copy
        echo "\\nSingle chunk detected, copying..." >> ${prefix}.merge.log
        cp \$(cat vcf_list.txt) ${prefix}.merged.vcf.gz
        cp \$(cat vcf_list.txt).tbi ${prefix}.merged.vcf.gz.tbi
    else
        # Multiple files, concatenate
        echo "\\nConcatenating chunks..." >> ${prefix}.merge.log
        bcftools concat \\
            --file-list vcf_list.txt \\
            --allow-overlaps \\
            --remove-duplicates \\
            --output-type z \\
            --output ${prefix}.merged.vcf.gz \\
            $args \\
            2>&1 | tee -a ${prefix}.merge.log
        
        # Index the output
        bcftools index -t ${prefix}.merged.vcf.gz
    fi
    
    # Get statistics
    echo "\\nMerge Statistics:" >> ${prefix}.merge.log
    n_samples=\$(bcftools query -l ${prefix}.merged.vcf.gz | wc -l)
    n_variants=\$(bcftools view -H ${prefix}.merged.vcf.gz | wc -l)
    echo "Total samples: \$n_samples" >> ${prefix}.merge.log
    echo "Total variants: \$n_variants" >> ${prefix}.merge.log
    
    # Check for chromosome consistency
    echo "\\nChromosomes in merged file:" >> ${prefix}.merge.log
    bcftools query -f '%CHROM\\n' ${prefix}.merged.vcf.gz | sort -u >> ${prefix}.merge.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi
    touch ${prefix}.merge.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}