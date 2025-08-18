process MERGE_FREQ_COMPARISON {
    tag "$meta.id"
    label 'process_low'
    
    container 'mamana/python-plotting:1.0.0'
    
    input:
    tuple val(meta), val(ref_name), path(imputed_freq), path(ref_freq)
    
    output:
    tuple val(meta), val(ref_name), path("*.freq_comparison.tsv"), emit: comparison
    path "*.summary.txt"                                          , emit: summary
    path "versions.yml"                                           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import numpy as np
    
    # Read imputed frequencies
    imp_df = pd.read_csv("${imputed_freq}", sep='\\t')
    
    # Read reference frequencies
    ref_df = pd.read_csv("${ref_freq}", sep='\\t')
    
    # Create a key for merging (CHROM:POS:REF:ALT)
    imp_df['key'] = imp_df['CHROM'].astype(str) + ':' + imp_df['POS'].astype(str) + ':' + \
                     imp_df['REF'] + ':' + imp_df['ALT']
    
    if not ref_df.empty:
        ref_df['key'] = ref_df['CHROM'].astype(str) + ':' + ref_df['POS'].astype(str) + ':' + \
                        ref_df['REF'] + ':' + ref_df['ALT']
        
        # Merge on the key
        merged_df = pd.merge(imp_df, ref_df[['key', 'REF_AF', 'REF_MAF']], 
                            on='key', how='left')
        
        # Fill missing reference values with NA
        merged_df['REF_AF'] = merged_df['REF_AF'].fillna('NA')
        merged_df['REF_MAF'] = merged_df['REF_MAF'].fillna('NA')
        
        # Calculate differences where both values are available
        mask = (merged_df['REF_AF'] != 'NA')
        merged_df['AF_DIFF'] = np.where(mask, 
                                        merged_df['IMP_AF'].astype(float) - merged_df['REF_AF'].astype(float),
                                        'NA')
        merged_df['MAF_DIFF'] = np.where(mask,
                                         merged_df['IMP_MAF'].astype(float) - merged_df['REF_MAF'].astype(float),
                                         'NA')
    else:
        # No reference data available
        merged_df = imp_df.copy()
        merged_df['REF_AF'] = 'NA'
        merged_df['REF_MAF'] = 'NA'
        merged_df['AF_DIFF'] = 'NA'
        merged_df['MAF_DIFF'] = 'NA'
    
    # Select and reorder columns for output
    output_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 
                   'IMP_AF', 'IMP_MAF', 'IMP_R2',
                   'REF_AF', 'REF_MAF', 
                   'AF_DIFF', 'MAF_DIFF']
    
    final_df = merged_df[output_cols]
    
    # Write comparison file
    output_file = "${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv"
    final_df.to_csv(output_file, sep='\\t', index=False)
    
    # Generate summary statistics
    summary_file = "${prefix}_${ref_name}_${chunk_id}.freq_comparison.summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Frequency Comparison Summary\\n")
        f.write(f"============================\\n")
        f.write(f"Total variants in imputed data: {len(imp_df)}\\n")
        
        if not ref_df.empty:
            f.write(f"Total variants in reference panel: {len(ref_df)}\\n")
            
            # Count matched variants
            matched = merged_df[merged_df['REF_AF'] != 'NA']
            f.write(f"Variants found in both: {len(matched)}\\n")
            f.write(f"Variants only in imputed: {len(merged_df) - len(matched)}\\n")
            
            if len(matched) > 0:
                # Calculate correlation if there are matched variants
                imp_af = matched['IMP_AF'].astype(float)
                ref_af = matched['REF_AF'].astype(float)
                correlation = np.corrcoef(imp_af, ref_af)[0, 1]
                
                # Calculate mean differences
                af_diff = matched['AF_DIFF'].astype(float)
                maf_diff = matched['MAF_DIFF'].astype(float)
                
                f.write(f"\\nFrequency Comparison Metrics:\\n")
                f.write(f"AF Correlation: {correlation:.4f}\\n")
                f.write(f"Mean AF difference: {af_diff.mean():.4f} (SD: {af_diff.std():.4f})\\n")
                f.write(f"Mean MAF difference: {maf_diff.mean():.4f} (SD: {maf_diff.std():.4f})\\n")
                f.write(f"Max AF difference: {af_diff.abs().max():.4f}\\n")
                f.write(f"Max MAF difference: {maf_diff.abs().max():.4f}\\n")
        else:
            f.write(f"No reference panel data available for comparison\\n")
    
    print(f"Merged frequency comparison for {len(final_df)} variants")
    
    # Create versions file
    with open('versions.yml', 'w') as f:
        f.write('${task.process}:\\n')
        f.write(f'    python: {pd.__version__}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_id = meta.chunk ?: ''
    """
    echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tIMP_AF\\tIMP_MAF\\tIMP_R2\\tREF_AF\\tREF_MAF\\tAF_DIFF\\tMAF_DIFF" > ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
    echo -e "chr21\\t10000001\\trs123\\tA\\tG\\t0.05\\t0.05\\t0.95\\t0.04\\t0.04\\t0.01\\t0.01" >> ${prefix}_${ref_name}_${chunk_id}.freq_comparison.tsv
    
    echo "Summary statistics" > ${prefix}_${ref_name}_${chunk_id}.freq_comparison.summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9
    END_VERSIONS
    """
}