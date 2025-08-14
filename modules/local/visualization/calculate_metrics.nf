process CALCULATE_IMPUTATION_METRICS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(imputed_vcf), path(info)
    tuple val(meta_orig), path(original_vcf)

    output:
    tuple val(meta), path("*.metrics.json"), emit: metrics
    tuple val(meta), path("*.concordance.tsv"), emit: concordance_data
    tuple val(meta), path("*.r2_data.tsv"), emit: r2_data
    tuple val(meta), path("*.summary.txt"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import json
    import gzip
    import numpy as np
    import pandas as pd
    from pathlib import Path
    
    # Calculate imputation metrics
    calculate_imputation_metrics.py \\
        --imputed-vcf ${imputed_vcf} \\
        --info-file ${info} \\
        --original-vcf ${original_vcf} \\
        --output-prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}
    
    # Generate summary statistics
    generate_summary_stats.py \\
        --metrics ${prefix}.metrics.json \\
        --output ${prefix}.summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"r2_mean": 0.95, "info_mean": 0.92}' > ${prefix}.metrics.json
    echo -e "MAF\\tConcordance\\n0.01\\t0.98" > ${prefix}.concordance.tsv
    echo -e "Position\\tR2\\n1000\\t0.95" > ${prefix}.r2_data.tsv
    echo "Summary: Mean R2=0.95" > ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        numpy: 1.24.0
        pandas: 2.0.0
    END_VERSIONS
    """
}