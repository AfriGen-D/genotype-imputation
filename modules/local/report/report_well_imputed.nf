process REPORT_WELL_IMPUTED {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(well_info), path(acc_info)
    
    output:
    tuple val(meta), val(ref_name), path("*.well_imputed.txt"), path("*.well_imputed_summary.txt"), emit: report
    path "versions.yml"                                                                            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import sys
    import numpy as np
    
    well_info_file = "${well_info}"
    report_file = "${prefix}_${ref_name}.well_imputed.txt"
    summary_file = "${prefix}_${ref_name}.well_imputed_summary.txt"
    
    # Process well imputed variants
    maf_bins = [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]
    bin_counts = {bin_range: 0 for bin_range in maf_bins}
    total_variants = 0
    
    with open(well_info_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\\t')
                if len(parts) >= 5:
                    try:
                        maf = float(parts[4])  # Assuming MAF is in column 5
                        total_variants += 1
                        
                        for bin_range in maf_bins:
                            if bin_range[0] <= maf < bin_range[1]:
                                bin_counts[bin_range] += 1
                                break
                    except (ValueError, IndexError):
                        continue
    
    # Write detailed report
    with open(report_file, 'w') as f:
        f.write("MAF_BIN\\tCOUNT\\tPERCENTAGE\\n")
        for bin_range, count in sorted(bin_counts.items()):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            f.write(f"{bin_range[0]:.2f}-{bin_range[1]:.2f}\\t{count}\\t{percentage:.2f}\\n")
    
    # Write summary
    with open(summary_file, 'w') as f:
        f.write(f"Total well imputed variants: {total_variants}\\n")
        f.write(f"Reference panel: ${ref_name}\\n")
        f.write(f"Sample: ${meta.id}\\n")
    
    print(f"Report generated: {report_file}")
    print(f"Summary generated: {summary_file}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.well_imputed.txt
    touch ${prefix}_${ref_name}.well_imputed_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}