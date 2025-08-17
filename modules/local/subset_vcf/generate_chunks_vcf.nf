process GENERATE_CHUNKS_VCF {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.chunks.txt"), emit: chunks
    path "*.summary.txt"                 , emit: summary
    path "versions.yml"                  , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_size = params.chunk_size ?: 5000000
    """
    #!/usr/bin/env python3
    
    import subprocess
    import sys
    
    # Get variant positions from VCF and determine ranges
    cmd = f"bcftools query -f '%CHROM\\t%POS\\n' ${vcf}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    chunks = []
    summary = []
    
    # Parse variant positions to get first and last per chromosome
    chrom_ranges = {}
    for line in result.stdout.strip().split('\\n'):
        if line:
            parts = line.split('\\t')
            if len(parts) >= 2:
                chrom = parts[0]
                pos = int(parts[1])
                if chrom not in chrom_ranges:
                    chrom_ranges[chrom] = [pos, pos]
                else:
                    chrom_ranges[chrom][0] = min(chrom_ranges[chrom][0], pos)
                    chrom_ranges[chrom][1] = max(chrom_ranges[chrom][1], pos)
    
    # Generate chunks based on actual variant range
    for chrom, (first_pos, last_pos) in chrom_ranges.items():
        chunk_count = 0
        # Generate chunks of specified size within the actual data range
        for start in range(first_pos, last_pos + 1, ${chunk_size}):
            end = min(start + ${chunk_size} - 1, last_pos)
            
            # Check if chunk contains any variants using bcftools
            cmd_check = f"bcftools view -H -r {chrom}:{start}-{end} ${vcf} | head -1 | wc -l"
            check_result = subprocess.run(cmd_check, shell=True, capture_output=True, text=True)
            has_variants = int(check_result.stdout.strip()) > 0
            
            if has_variants:
                chunk_id = f"chunk_{chrom}_{start}_{end}"
                chunks.append(f"{chrom}\\t{start}\\t{end}\\t{chunk_id}")
                chunk_count += 1
        
        summary.append(f"{chrom}: {chunk_count} chunks (range: {first_pos}-{last_pos})")
    
    # Write chunks file
    with open("${prefix}.chunks.txt", "w") as f:
        for chunk in chunks:
            f.write(chunk + "\\n")
    
    # Write summary file
    with open("${prefix}.summary.txt", "w") as f:
        f.write(f"Total chunks: {len(chunks)}\\n")
        for line in summary:
            f.write(line + "\\n")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    python: ' + sys.version.split()[0] + '\\n')
        
    import subprocess
    bcftools_version = subprocess.run(['bcftools', '--version'], capture_output=True, text=True).stdout.split('\\n')[0].split()[-1]
    with open("versions.yml", "a") as f:
        f.write(f'    bcftools: {bcftools_version}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chunks.txt
    touch ${prefix}.summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        bcftools: 1.20
    END_VERSIONS
    """
}