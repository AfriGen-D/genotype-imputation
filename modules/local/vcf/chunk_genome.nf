process CHUNK_GENOME {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)
    val chunk_size

    output:
    tuple val(meta), path("${prefix}_chunks.txt"), emit: chunk_file
    path("${prefix}_*.bed"), emit: chunk_beds
    val(meta), emit: chunks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import subprocess
    import pandas as pd
    from pathlib import Path
    
    # Get chromosome lengths from VCF
    cmd = f"bcftools index -s {vcf}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    chunks = []
    chunk_num = 0
    
    for line in result.stdout.strip().split('\\n'):
        parts = line.split('\\t')
        chrom = parts[0]
        length = int(parts[1])
        
        # Create chunks for this chromosome
        for start in range(1, length, ${chunk_size}):
            end = min(start + ${chunk_size} - 1, length)
            chunk_num += 1
            
            # Create BED file for chunk
            bed_file = f"${prefix}_chunk_{chunk_num:04d}.bed"
            with open(bed_file, 'w') as f:
                f.write(f"{chrom}\\t{start}\\t{end}\\n")
            
            chunks.append({
                'chunk_id': chunk_num,
                'chr': chrom,
                'start': start,
                'end': end,
                'region': f"{chrom}:{start}-{end}",
                'bed': bed_file
            })
    
    # Write chunks file
    chunks_df = pd.DataFrame(chunks)
    chunks_df.to_csv("${prefix}_chunks.txt", sep='\\t', index=False)
    
    print(f"Created {len(chunks)} chunks of size {${chunk_size}} bp")
    
    # Write versions file
    import sys
    with open("versions.yml", "w") as f:
        f.write(f'''${task.process}:
        python: {sys.version.split()[0]}
        pandas: {pd.__version__}
    ''')
    """
}