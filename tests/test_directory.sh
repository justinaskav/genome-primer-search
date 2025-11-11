#!/bin/bash
# Test script: Directory input with multiple genome files
# Tests the pipeline with a directory containing multiple .txt files
# Verifies: file merging, deduplication, parameter passing

set -e

# Parse command line arguments
ENABLE_THERMO=false
if [[ "$1" == "--with-thermo" ]]; then
    ENABLE_THERMO=true
fi

echo "========================================"
echo "Test: Directory Input (Multiple Files)"
echo "========================================"
echo "Testing with: tests/data/test_genomes_dir/"
echo ""
echo "Directory contains:"
ls -1 tests/data/test_genomes_dir/*.txt | sed 's/^/  - /'
echo ""
echo "Expected behavior:"
echo "  ✓ Merge all .txt files"
echo "  ✓ Deduplicate genome entries"
echo "  ✓ Filter empty lines and comments"
echo "  ✓ Process unique genomes: 4 total"
echo "    (GCA_000005845.2, Escherichia coli, GCA_000027305.1, GCA_000006885.1)"
echo ""

if [ "$ENABLE_THERMO" = true ]; then
    echo "With thermodynamic analysis enabled"
    echo ""
fi

# Clean previous test output
rm -rf test_output_directory

# Build Nextflow command
NF_CMD="nextflow run main.nf \
    --genomes tests/data/test_genomes_dir \
    --primers tests/data/test_primers.txt \
    --outdir test_output_directory \
    --min_amplicon_size 50 \
    --max_amplicon_size 2000 \
    -resume"

# Add thermodynamic analysis if requested
if [ "$ENABLE_THERMO" = true ]; then
    NF_CMD="$NF_CMD \
    --enable_thermo_analysis \
    --thermo_references tests/data/test_references.fasta \
    --thermo_master_mix DreamTaq \
    --thermo_annealing_temp 60"
fi

# Run pipeline with test data
echo "Running pipeline..."
eval $NF_CMD

echo ""
echo "========================================"
echo "✓ Directory test completed!"
echo "========================================"
echo ""
echo "Check results in: test_output_directory/"
echo "  - test_output_directory/results/genomes/      Downloaded genomes"
echo "  - test_output_directory/results/filtered/     Filtered amplicons"
echo "  - test_output_directory/reports/summary.html  HTML report (with parameters)"
echo "  - test_output_directory/reports/summary.tsv   TSV summary"

if [ "$ENABLE_THERMO" = true ]; then
    echo "  - test_output_directory/thermo_reports/       Thermodynamic analysis"
fi

echo ""
echo "Verification steps:"
echo "  1. Check genome count:"
GENOME_COUNT=$(find test_output_directory/results/genomes -name "*.fna" 2>/dev/null | wc -l | tr -d ' ')
echo "     → Genomes downloaded: $GENOME_COUNT (expected: 4)"

echo ""
echo "  2. Check summary.html for:"
echo "     → Filter criteria: '50 - 2,000 bp'"
echo "     → Total genomes: 4 (deduplicated from 6 entries)"
echo ""
echo "  3. Verify no duplicates:"
echo "     → Count genome files in results/genomes/"
echo "     → Should have exactly 4 unique genomes"
echo ""
