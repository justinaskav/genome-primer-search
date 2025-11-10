#!/bin/bash
# Test script for Genome Primer Search Pipeline
# Runs a quick validation test with small genomes

set -e

echo "================================"
echo "Running pipeline test..."
echo "================================"
echo ""

# Run pipeline with test data
nextflow run main.nf \
    --genomes tests/data/test_genomes.txt \
    --primers tests/data/test_primers.txt \
    --outdir test_output \
    -resume

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
