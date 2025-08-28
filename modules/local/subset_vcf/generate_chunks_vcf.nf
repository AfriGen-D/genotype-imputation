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
    def chunk_size = params.chunk_size ?: 50000000  // Default chunk size
    def min_variants = params.min_chunk_variants ?: 1000  // Minimum variants per chunk
    def max_chunk_size = params.max_chunk_size ?: 100000000  // Maximum chunk size (100Mb)
    """
    #!/usr/bin/env python3
    
    import subprocess
    import sys
    
    # Configuration
    BASE_CHUNK_SIZE = ${chunk_size}  # Base chunk size in bp
    MIN_VARIANTS = ${min_variants}    # Minimum number of variants required per chunk
    MAX_CHUNK_SIZE = ${max_chunk_size}  # Maximum chunk size to prevent memory issues
    
    print(f"Adaptive chunking configuration:")
    print(f"  Base chunk size: {BASE_CHUNK_SIZE:,} bp")
    print(f"  Min variants per chunk: {MIN_VARIANTS:,}")
    print(f"  Max chunk size: {MAX_CHUNK_SIZE:,} bp")
    print()
    
    # Get all variant positions from VCF
    cmd = f"bcftools query -f '%CHROM\\t%POS\\n' ${vcf}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    chunks = []
    summary = []
    
    # Parse variant positions per chromosome
    chrom_positions = {}
    for line in result.stdout.strip().split('\\n'):
        if line:
            parts = line.split('\\t')
            if len(parts) >= 2:
                chrom = parts[0]
                pos = int(parts[1])
                if chrom not in chrom_positions:
                    chrom_positions[chrom] = []
                chrom_positions[chrom].append(pos)
    
    # Generate adaptive chunks for each chromosome
    for chrom, positions in chrom_positions.items():
        if not positions:
            continue
            
        positions.sort()
        first_pos = positions[0]
        last_pos = positions[-1]
        total_variants = len(positions)
        
        print(f"Processing {chrom}: {total_variants:,} variants in range {first_pos:,}-{last_pos:,}")
        
        chunk_count = 0
        current_start = first_pos
        
        while current_start <= last_pos:
            # Start with base chunk size
            current_end = min(current_start + BASE_CHUNK_SIZE - 1, last_pos)
            
            # Count variants in proposed chunk
            variants_in_chunk = sum(1 for p in positions if current_start <= p <= current_end)
            
            # Expand chunk if it has too few variants (but respect max size)
            while variants_in_chunk < MIN_VARIANTS and current_end < last_pos:
                # Extend chunk by 25% or to next significant variant cluster
                extension = min(BASE_CHUNK_SIZE // 4, MAX_CHUNK_SIZE - (current_end - current_start))
                new_end = min(current_end + extension, last_pos)
                
                # Check if we'd exceed max chunk size
                if (new_end - current_start) > MAX_CHUNK_SIZE:
                    break
                    
                current_end = new_end
                variants_in_chunk = sum(1 for p in positions if current_start <= p <= current_end)
                
                # If still not enough variants and at end of chromosome, accept what we have
                if current_end >= last_pos:
                    break
            
            # Only create chunk if it has any variants
            if variants_in_chunk > 0:
                chunk_id = f"chunk_{chrom}_{current_start}_{current_end}"
                chunks.append(f"{chrom}\\t{current_start}\\t{current_end}\\t{chunk_id}")
                chunk_count += 1
                chunk_size_mb = (current_end - current_start + 1) / 1000000
                print(f"  Chunk {chunk_count}: {chrom}:{current_start:,}-{current_end:,} ({chunk_size_mb:.1f} Mb, {variants_in_chunk:,} variants)")
                
                # If chunk was expanded due to low variant density, warn
                if chunk_size_mb > (BASE_CHUNK_SIZE / 1000000 * 1.5):
                    print(f"    ⚠ Expanded chunk size due to low variant density")
            
            # Move to next chunk
            current_start = current_end + 1
        
        summary.append(f"{chrom}: {chunk_count} chunks (range: {first_pos:,}-{last_pos:,}, {total_variants:,} variants)")
    
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