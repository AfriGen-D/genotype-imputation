process AGGREGATE_CHUNKS_TO_CHR {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/chromosome/${meta.dataset}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(chunk_reports), path(chunk_info)
    
    output:
    tuple val(meta), val(ref_name), path("*.chr_summary.json"), emit: chr_summary
    tuple val(meta), val(ref_name), path("*.chr_stats.txt")   , emit: chr_stats
    path "versions.yml"                                        , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Aggregate chunk-level metrics to chromosome level
    aggregate_to_chromosome.py \\
        --chunk-reports ${chunk_reports} \\
        --chunk-info ${chunk_info} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        --chromosome ${meta.chromosome} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.chr_summary.json
    touch ${prefix}_${ref_name}.chr_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}