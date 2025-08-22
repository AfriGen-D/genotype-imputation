process MERGE_ADJACENT_CHUNKS {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/vcf-processing:bcftools-1.20'
    
    input:
    tuple val(meta), val(adjacent_meta), path(main_vcf), path(adjacent_vcf)
    
    output:
    tuple val(meta), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: vcf
    path "*.merge.log"                                                    , emit: log
    path "versions.yml"                                                  , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Merging failed chunk variants into adjacent passing chunk" > ${prefix}.merge.log
    echo "Main chunk: ${meta.id} (${meta.contig}:${meta.start}-${meta.end})" >> ${prefix}.merge.log
    echo "Failed adjacent chunk: ${adjacent_meta.id} (${adjacent_meta.contig}:${adjacent_meta.start}-${adjacent_meta.end})" >> ${prefix}.merge.log
    
    # Ensure BCF files have indexes (create if missing)
    if [ ! -f "${main_vcf}.csi" ]; then
        echo "Creating index for main chunk BCF" >> ${prefix}.merge.log
        bcftools index ${main_vcf}
    fi
    if [ ! -f "${adjacent_vcf}.csi" ]; then
        echo "Creating index for adjacent chunk BCF" >> ${prefix}.merge.log
        bcftools index ${adjacent_vcf}
    fi
    
    # Count variants before merge
    main_count=\$(bcftools view -H ${main_vcf} | wc -l)
    adjacent_count=\$(bcftools view -H ${adjacent_vcf} | wc -l)
    echo "Variants in main chunk: \$main_count" >> ${prefix}.merge.log
    echo "Variants in failed chunk: \$adjacent_count" >> ${prefix}.merge.log
    
    # Merge the chunks, removing duplicates
    bcftools concat \\
        ${main_vcf} \\
        ${adjacent_vcf} \\
        --allow-overlaps \\
        --remove-duplicates \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        --threads ${task.cpus}
    
    # Index the output
    bcftools index -t ${prefix}.merged.vcf.gz
    
    # Report merge statistics
    merged_count=\$(bcftools view -H ${prefix}.merged.vcf.gz | wc -l)
    echo "Variants after merge: \$merged_count" >> ${prefix}.merge.log
    echo "Variants added from failed chunk: \$((merged_count - main_count))" >> ${prefix}.merge.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi
    touch ${prefix}.merge.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}