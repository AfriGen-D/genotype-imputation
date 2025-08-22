#!/bin/bash
# Regenerate JSON summaries and plots for existing results

set -e

RESULTS_DIR="/scratch3/users/mamana/results/reports"
SCRIPT_DIR="/users/mamana/genotype-imputation/bin"

echo "Step 1: Generating missing chunk JSON files..."

# Find all chunk directories that have reports but no JSON
for chunk_dir in $RESULTS_DIR/awigen_500_b38_chr*; do
    if [ -d "$chunk_dir" ]; then
        chunk_name=$(basename $chunk_dir)
        
        # Check if JSON already exists
        if [ ! -f "$chunk_dir/${chunk_name}_H3AR6x.chunk_summary.json" ]; then
            # Check if required files exist
            if [ -f "$chunk_dir/${chunk_name}_H3AR6x.accuracy.txt" ] && \
               [ -f "$chunk_dir/${chunk_name}_H3AR6x.well_imputed.txt" ] && \
               [ -f "$chunk_dir/${chunk_name}_H3AR6x.well_imputed_summary.txt" ]; then
                
                echo "  Generating JSON for $chunk_name..."
                cd "$chunk_dir"
                python3 "$SCRIPT_DIR/generate_chunk_json.py" \
                    --accuracy-file "${chunk_name}_H3AR6x.accuracy.txt" \
                    --well-imputed-file "${chunk_name}_H3AR6x.well_imputed.txt" \
                    --summary-file "${chunk_name}_H3AR6x.well_imputed_summary.txt" \
                    --chunk-id "$chunk_name" \
                    --ref-name "H3AR6x" \
                    --output "${chunk_name}_H3AR6x.chunk_summary.json" 2>/dev/null || true
            fi
        fi
    fi
done

echo "Step 2: Re-aggregating chromosome-level summaries..."

# For each chromosome, aggregate chunk JSONs
for chr in {1..22}; do
    echo "  Processing chromosome $chr..."
    
    # Find all chunk JSONs for this chromosome
    chunk_jsons=$(find $RESULTS_DIR -name "awigen_500_b38_chr${chr}_*_H3AR6x.chunk_summary.json" -type f 2>/dev/null)
    
    if [ ! -z "$chunk_jsons" ]; then
        # Create chromosome directory if it doesn't exist
        chr_dir="$RESULTS_DIR/chromosome/awigen_500_b38"
        mkdir -p "$chr_dir"
        
        cd "$chr_dir"
        python3 "$SCRIPT_DIR/aggregate_to_chromosome.py" \
            --chunk-reports $chunk_jsons \
            --output-prefix "awigen_500_b38_chr${chr}" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "chr${chr}" 2>/dev/null || true
    fi
done

echo "Step 3: Re-generating chromosome plots..."

# Generate plots for each chromosome summary
chr_dir="$RESULTS_DIR/chromosome/awigen_500_b38"
plot_dir="$chr_dir/plots"
mkdir -p "$plot_dir"

for chr_json in $chr_dir/*_H3AR6x.chr_summary.json; do
    if [ -f "$chr_json" ]; then
        chr_name=$(basename "$chr_json" .chr_summary.json | sed 's/_H3AR6x//')
        chr_num=$(echo "$chr_name" | grep -oE 'chr[0-9]+')
        
        echo "  Generating plots for $chr_num..."
        
        cd "$plot_dir"
        
        # Generate accuracy plot
        python3 "$SCRIPT_DIR/plot_chr_accuracy.py" \
            --chr-summary "$chr_json" \
            --output-prefix "$chr_name" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "$chr_num" 2>/dev/null || true
        
        # Generate performance plot
        python3 "$SCRIPT_DIR/plot_chr_performance.py" \
            --chr-summary "$chr_json" \
            --output-prefix "$chr_name" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "$chr_num" 2>/dev/null || true
        
        # Generate MAF analysis plot
        python3 "$SCRIPT_DIR/plot_chr_maf_analysis.py" \
            --chr-summary "$chr_json" \
            --output-prefix "$chr_name" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "$chr_num" 2>/dev/null || true
    fi
done

echo "Step 4: Aggregating genome-level summary..."

# Aggregate all chromosome summaries to genome level
genome_dir="$RESULTS_DIR/genome/awigen_500_b38"
mkdir -p "$genome_dir"

chr_jsons=$(find $chr_dir -name "*_H3AR6x.chr_summary.json" -type f)

if [ ! -z "$chr_jsons" ]; then
    cd "$genome_dir"
    python3 "$SCRIPT_DIR/aggregate_to_genome.py" \
        --chr-summaries $chr_jsons \
        --output-prefix "awigen_500_b38" \
        --ref-name "H3AR6x" \
        --dataset "awigen_500_b38" 2>/dev/null || true
fi

echo "Step 5: Generating genome-level plots..."

genome_json="$genome_dir/awigen_500_b38_H3AR6x.genome_summary.json"

if [ -f "$genome_json" ]; then
    cd "$genome_dir"
    
    python3 "$SCRIPT_DIR/plot_genome_summary.py" \
        --genome-summary "$genome_json" \
        --output-prefix "awigen_500_b38" \
        --ref-name "H3AR6x" \
        --dataset "awigen_500_b38" 2>/dev/null || true
fi

echo "Report regeneration complete!"
echo "Generated files:"
echo "  - Chunk JSONs: $(find $RESULTS_DIR -name "*chunk_summary.json" | wc -l)"
echo "  - Chromosome JSONs: $(find $chr_dir -name "*.chr_summary.json" | wc -l)"
echo "  - Chromosome plots: $(find $plot_dir -name "*.pdf" | wc -l)"