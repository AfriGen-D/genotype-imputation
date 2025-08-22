#!/bin/bash
# Regenerate plots using the plotting container

set -e

RESULTS_DIR="/scratch3/users/mamana/results/reports"
SCRIPT_DIR="/users/mamana/genotype-imputation/bin"
CONTAINER="mamana/python-plotting:1.0.0"

echo "Regenerating plots using container: $CONTAINER"

# Use singularity if available, otherwise docker
if command -v singularity &> /dev/null; then
    RUNNER="singularity exec docker://$CONTAINER"
elif command -v docker &> /dev/null; then
    RUNNER="docker run --rm -v $RESULTS_DIR:$RESULTS_DIR -v $SCRIPT_DIR:$SCRIPT_DIR -w /tmp $CONTAINER"
else
    echo "Error: Neither singularity nor docker found!"
    exit 1
fi

echo "Using: $RUNNER"

echo "Regenerating chromosome plots..."

# Generate plots for each chromosome summary
chr_dir="$RESULTS_DIR/chromosome/awigen_500_b38"
plot_dir="$chr_dir/plots"
mkdir -p "$plot_dir"

for chr in {1..22}; do
    chr_json="$chr_dir/awigen_500_b38_chr${chr}_H3AR6x.chr_summary.json"
    if [ -f "$chr_json" ]; then
        echo "  Generating plots for chr${chr}..."
        
        # Accuracy plot
        $RUNNER python3 "$SCRIPT_DIR/plot_chr_accuracy.py" \
            --chr-summary "$chr_json" \
            --output-prefix "$plot_dir/awigen_500_b38_chr${chr}" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "chr${chr}" || echo "    Failed: accuracy plot"
        
        # Performance plot
        $RUNNER python3 "$SCRIPT_DIR/plot_chr_performance.py" \
            --chr-summary "$chr_json" \
            --output-prefix "$plot_dir/awigen_500_b38_chr${chr}" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "chr${chr}" || echo "    Failed: performance plot"
        
        # MAF analysis plot
        $RUNNER python3 "$SCRIPT_DIR/plot_chr_maf_analysis.py" \
            --chr-summary "$chr_json" \
            --output-prefix "$plot_dir/awigen_500_b38_chr${chr}" \
            --ref-name "H3AR6x" \
            --dataset "awigen_500_b38" \
            --chromosome "chr${chr}" || echo "    Failed: MAF analysis plot"
    fi
done

echo "Regenerating genome-level plots..."

genome_dir="$RESULTS_DIR/genome/awigen_500_b38"
genome_json="$genome_dir/awigen_500_b38_H3AR6x.genome_summary.json"

if [ -f "$genome_json" ]; then
    $RUNNER python3 "$SCRIPT_DIR/plot_genome_summary.py" \
        --genome-summary "$genome_json" \
        --output-prefix "$genome_dir/awigen_500_b38" \
        --ref-name "H3AR6x" \
        --dataset "awigen_500_b38" || echo "Failed: genome summary plot"
    
    $RUNNER python3 "$SCRIPT_DIR/plot_genome_accuracy.py" \
        --genome-summary "$genome_json" \
        --output-prefix "$genome_dir/awigen_500_b38" \
        --ref-name "H3AR6x" \
        --dataset "awigen_500_b38" || echo "Failed: genome accuracy plot"
    
    $RUNNER python3 "$SCRIPT_DIR/plot_genome_performance.py" \
        --genome-summary "$genome_json" \
        --output-prefix "$genome_dir/awigen_500_b38" \
        --ref-name "H3AR6x" \
        --dataset "awigen_500_b38" || echo "Failed: genome performance plot"
fi

echo "Plot regeneration complete!"
echo "Chromosome plots: $(find $plot_dir -name "*.pdf" -newer /tmp/test_chr1_H3AR6x.chr_summary.json 2>/dev/null | wc -l) new"
echo "Genome plots: $(find $genome_dir -name "*.pdf" -newer /tmp/test_chr1_H3AR6x.chr_summary.json 2>/dev/null | wc -l) new"