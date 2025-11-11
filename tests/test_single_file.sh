#!/bin/bash
# Test script: Single genome file input
# Tests the pipeline with a single .txt file containing genome names/accessions

set -e

# Parse command line arguments
ENABLE_THERMO=false
if [[ "$1" == "--with-thermo" ]]; then
    ENABLE_THERMO=true
fi

echo "========================================"
echo "Test: Single Genome File Input"
echo "========================================"
echo "Testing with: tests/data/test_genomes.txt"
if [ "$ENABLE_THERMO" = true ]; then
    echo "With thermodynamic analysis enabled"
fi
echo ""

# Clean previous test output
rm -rf test_output_single

# Build Nextflow command
NF_CMD="nextflow run main.nf \
    --genomes tests/data/test_genomes.txt \
    --primers tests/data/test_primers.txt \
    --outdir test_output_single \
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
echo "✓ Single file test completed!"
echo "========================================"
echo ""
echo "Check results in: test_output_single/"
echo "  - test_output_single/results/genomes/      Downloaded genomes"
echo "  - test_output_single/results/filtered/     Filtered amplicons"
echo "  - test_output_single/reports/summary.html  HTML report (with parameters)"
echo "  - test_output_single/reports/summary.tsv   TSV summary"

if [ "$ENABLE_THERMO" = true ]; then
    echo "  - test_output_single/thermo_reports/       Thermodynamic analysis"
fi

echo ""
echo "Verification:"
echo "  → Check that summary.html shows: '50 - 2,000 bp' filter criteria"
echo "  → Verify no duplicate genomes were processed"
echo ""
