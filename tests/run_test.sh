#!/bin/bash
# Test script for Genome Primer Search Pipeline
# Runs a quick validation test with small genomes

set -e

# Parse command line arguments
ENABLE_THERMO=false
if [[ "$1" == "--with-thermo" ]]; then
    ENABLE_THERMO=true
fi

echo "================================"
echo "Running pipeline test..."
if [ "$ENABLE_THERMO" = true ]; then
    echo "With thermodynamic analysis"
fi
echo "================================"
echo ""

# Build Nextflow command
NF_CMD="nextflow run main.nf \
    --genomes tests/data/test_genomes.txt \
    --primers tests/data/test_primers.txt \
    --outdir test_output \
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
eval $NF_CMD

echo ""
echo "================================"
echo "Test completed successfully!"
echo "================================"
echo ""
echo "Check results in: test_output/"
echo "  - test_output/results/genomes/      Downloaded genomes"
echo "  - test_output/results/filtered/     Filtered amplicons"
echo "  - test_output/reports/summary.html  HTML report"
echo "  - test_output/reports/summary.tsv   TSV summary"

if [ "$ENABLE_THERMO" = true ]; then
    echo "  - test_output/thermo_reports/       Thermodynamic analysis"
fi

echo ""
echo "Usage:"
echo "  bash tests/run_test.sh              # Basic test"
echo "  bash tests/run_test.sh --with-thermo  # Test with thermodynamic analysis"
