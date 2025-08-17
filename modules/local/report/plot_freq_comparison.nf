process PLOT_FREQ_COMPARISON {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(imputed_vcf), path(vcf_index)
    
    output:
    tuple val(meta), val(ref_name), path("*.freq_comparison.png"), emit: plot
    path "versions.yml"                                          , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    import subprocess
    import sys
    
    # Extract allele frequencies using bcftools
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' ${imputed_vcf}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Parse frequencies
    frequencies = []
    positions = []
    
    for line in result.stdout.strip().split('\\n'):
        if line:
            parts = line.split('\\t')
            if len(parts) >= 5:
                try:
                    pos = int(parts[1])
                    af = float(parts[4]) if parts[4] != '.' else 0.5
                    positions.append(pos)
                    frequencies.append(af)
                except (ValueError, IndexError):
                    continue
    
    # Create frequency distribution plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Histogram of allele frequencies
    if frequencies:
        ax1.hist(frequencies, bins=50, edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Allele Frequency')
        ax1.set_ylabel('Count')
        ax1.set_title('Distribution of Imputed Allele Frequencies')
        ax1.grid(True, alpha=0.3)
        
        # Frequency vs position (if we have position info)
        if positions and len(positions) == len(frequencies):
            ax2.scatter(positions[:1000], frequencies[:1000], alpha=0.5, s=5)  # Limit to first 1000 for visibility
            ax2.set_xlabel('Genomic Position')
            ax2.set_ylabel('Allele Frequency')
            ax2.set_title('Allele Frequency Along Chromosome')
            ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f'Frequency Analysis - ${meta.id} (${ref_name})', fontsize=14)
    plt.tight_layout()
    plt.savefig("${prefix}_${ref_name}.freq_comparison.png", dpi=150)
    plt.close()
    
    print(f"Frequency comparison plot saved: ${prefix}_${ref_name}.freq_comparison.png")
    print(f"Total variants analyzed: {len(frequencies)}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${ref_name}.freq_comparison.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}