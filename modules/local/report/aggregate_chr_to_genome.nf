process AGGREGATE_CHR_TO_GENOME {
    tag "$meta.id"
    label 'process_medium'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/genome/${meta.dataset}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(chr_summaries)
    
    output:
    tuple val(meta), val(ref_name), path("*.genome_summary.json"), emit: genome_summary
    tuple val(meta), val(ref_name), path("*.genome_stats.txt")   , emit: genome_stats
    path "versions.yml"                                           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Aggregate chromosome-level summaries to genome-wide level
    aggregate_to_genome.py \\
        --chr-summaries ${chr_summaries} \\
        --output-prefix ${prefix} \\
        --ref-name ${ref_name} \\
        --dataset ${meta.dataset} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.genome_summary.json
    touch ${prefix}_${ref_name}.genome_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}