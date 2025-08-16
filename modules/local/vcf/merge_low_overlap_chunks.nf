process MERGE_LOW_OVERLAP_CHUNKS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(chunk_vcfs), path(chunk_indices), path(overlap_matrix)
    val min_overlap_ratio

    output:
    tuple val(meta), path("*_merged_chunk_*.vcf.gz"), path("*_merged_chunk_*.vcf.gz.tbi"), emit: merged_chunks
    path "chunk_merge_report.txt", emit: merge_report
    path "merged_chunk_list.txt", emit: chunk_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def min_ratio_pct = min_overlap_ratio * 100  // Convert to percentage (0.00001 -> 0.001%)
    """
    echo "Analyzing chunks for merging based on minimum overlap ratio: ${min_ratio_pct}%"
    
    # Parse overlap matrix to identify chunks needing merging
    # Skip header and extract chunk names and overlap percentages
    tail -n +2 ${overlap_matrix} | while IFS=',' read -r chunk ref_panel target_vars ref_vars exact exact_pct pos pos_pct; do
        echo "\$chunk:\$exact_pct" >> chunk_overlaps.txt
    done
    
    # Sort chunks by name to maintain order
    sort -t: -k1,1 chunk_overlaps.txt > sorted_overlaps.txt
    
    # Create merge plan
    python3 <<'EOF'
import os
import subprocess

min_ratio = ${min_overlap_ratio} * 100  # Convert to percentage (0.00001 = 0.001%)
chunks_to_merge = []
current_group = []
chunk_files = {}

# Read chunk files and their overlap percentages
with open("sorted_overlaps.txt", "r") as f:
    for line in f:
        parts = line.strip().split(":")
        if len(parts) == 2:
            chunk_name = parts[0]
            overlap_pct = float(parts[1])
            
            # Find corresponding VCF file
            vcf_file = f"{chunk_name}.vcf.gz"
            if os.path.exists(vcf_file):
                chunk_files[chunk_name] = {
                    'file': vcf_file,
                    'overlap': overlap_pct
                }

# Sort chunks by their numeric ID
sorted_chunks = sorted(chunk_files.keys(), 
                      key=lambda x: int(x.split('_chunk_')[1]) if '_chunk_' in x else 0)

# Group chunks for merging
merge_groups = []
current_group = []

for i, chunk in enumerate(sorted_chunks):
    overlap = chunk_files[chunk]['overlap']
    
    if overlap < min_ratio:
        # This chunk needs to be merged
        if not current_group:
            # Start a new merge group
            # Include previous chunk if available and not already in a group
            if i > 0 and sorted_chunks[i-1] not in [c for g in merge_groups for c in g]:
                current_group.append(sorted_chunks[i-1])
        current_group.append(chunk)
        
        # Look ahead to see if next chunk also needs merging
        if i < len(sorted_chunks) - 1:
            next_overlap = chunk_files[sorted_chunks[i+1]]['overlap']
            if next_overlap >= min_ratio:
                # Next chunk has good overlap, finalize this group
                merge_groups.append(current_group)
                current_group = []
    else:
        # This chunk has good overlap
        if current_group:
            # Finalize the merge group
            merge_groups.append(current_group)
            current_group = []

# Add any remaining group
if current_group:
    merge_groups.append(current_group)

# Write merge plan
with open("merge_plan.txt", "w") as f:
    for i, group in enumerate(merge_groups):
        f.write(f"Group {i}: {','.join(group)}\\n")

# Execute merges
merged_chunks = []
chunk_id = 0

# First, add chunks that don't need merging
for chunk in sorted_chunks:
    if not any(chunk in group for group in merge_groups):
        # This chunk doesn't need merging, just copy it
        new_name = f"${prefix}_merged_chunk_{chunk_id:04d}"
        subprocess.run(f"cp {chunk_files[chunk]['file']} {new_name}.vcf.gz", shell=True)
        subprocess.run(f"cp {chunk_files[chunk]['file']}.tbi {new_name}.vcf.gz.tbi", shell=True)
        merged_chunks.append(new_name)
        chunk_id += 1

# Now merge the groups
for group_id, group in enumerate(merge_groups):
    if len(group) > 1:
        new_name = f"${prefix}_merged_chunk_{chunk_id:04d}"
        vcf_list = " ".join([chunk_files[c]['file'] for c in group])
        
        print(f"Merging {len(group)} chunks into {new_name}: {', '.join(group)}")
        
        # Concatenate VCFs in the group
        cmd = f"bcftools concat {vcf_list} -a -D -Oz -o {new_name}.vcf.gz"
        subprocess.run(cmd, shell=True, check=True)
        
        # Index the merged chunk
        subprocess.run(f"bcftools index -t {new_name}.vcf.gz", shell=True, check=True)
        
        merged_chunks.append(new_name)
        chunk_id += 1
    elif len(group) == 1:
        # Single chunk in group (shouldn't happen but handle it)
        new_name = f"${prefix}_merged_chunk_{chunk_id:04d}"
        subprocess.run(f"cp {chunk_files[group[0]]['file']} {new_name}.vcf.gz", shell=True)
        subprocess.run(f"cp {chunk_files[group[0]]['file']}.tbi {new_name}.vcf.gz.tbi", shell=True)
        merged_chunks.append(new_name)
        chunk_id += 1

# Write merged chunk list
with open("merged_chunk_list.txt", "w") as f:
    for chunk in sorted(merged_chunks):
        f.write(f"{chunk}\\n")

print(f"Created {len(merged_chunks)} merged chunks from {len(sorted_chunks)} original chunks")
EOF

    # Create merge report
    cat > chunk_merge_report.txt <<EOF
====================================
CHUNK MERGE REPORT
====================================
Minimum overlap ratio: ${min_ratio_pct}%
Original chunks: \$(ls *_chunk_*.vcf.gz 2>/dev/null | wc -l)
Merged chunks: \$(ls *_merged_chunk_*.vcf.gz 2>/dev/null | wc -l)

MERGE PLAN:
EOF
    
    if [ -f merge_plan.txt ]; then
        cat merge_plan.txt >> chunk_merge_report.txt
    fi
    
    echo "" >> chunk_merge_report.txt
    echo "MERGED CHUNKS:" >> chunk_merge_report.txt
    
    # Report statistics for each merged chunk
    for merged in *_merged_chunk_*.vcf.gz; do
        if [ -f "\$merged" ]; then
            n_vars=\$(bcftools index -n \$merged)
            echo "  \$merged: \$n_vars variants" >> chunk_merge_report.txt
        fi
    done
    
    # Display summary
    echo "Merge complete:"
    cat chunk_merge_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}