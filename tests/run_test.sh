#!/bin/bash
# Comprehensive Test Runner for Genome Primer Search Pipeline
# Runs multiple test scenarios to validate pipeline functionality

set -e

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse command line arguments
TEST_MODE="all"
ENABLE_THERMO=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --single)
            TEST_MODE="single"
            shift
            ;;
        --directory)
            TEST_MODE="directory"
            shift
            ;;
        --all)
            TEST_MODE="all"
            shift
            ;;
        --with-thermo)
            ENABLE_THERMO=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--single|--directory|--all] [--with-thermo]"
            exit 1
            ;;
    esac
done

echo ""
echo "=================================================="
echo "  Genome Primer Search Pipeline - Test Suite"
echo "=================================================="
echo ""

# Build thermo flag
THERMO_FLAG=""
if [ "$ENABLE_THERMO" = true ]; then
    THERMO_FLAG="--with-thermo"
    echo -e "${YELLOW}Thermodynamic analysis: ENABLED${NC}"
else
    echo "Thermodynamic analysis: disabled"
fi
echo ""

# Run tests based on mode
if [ "$TEST_MODE" = "all" ] || [ "$TEST_MODE" = "single" ]; then
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}TEST 1: Single Genome File Input${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    bash tests/test_single_file.sh $THERMO_FLAG

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Single file test PASSED${NC}"
    else
        echo "✗ Single file test FAILED"
        exit 1
    fi
    echo ""
fi

if [ "$TEST_MODE" = "all" ] || [ "$TEST_MODE" = "directory" ]; then
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}TEST 2: Directory Input (Multiple Files)${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    bash tests/test_directory.sh $THERMO_FLAG

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Directory test PASSED${NC}"
    else
        echo "✗ Directory test FAILED"
        exit 1
    fi
    echo ""
fi

echo ""
echo "=================================================="
echo -e "${GREEN}  ✓ All tests completed successfully!${NC}"
echo "=================================================="
echo ""
echo "Test outputs:"
if [ "$TEST_MODE" = "all" ] || [ "$TEST_MODE" = "single" ]; then
    echo "  - test_output_single/       Single file test results"
fi
if [ "$TEST_MODE" = "all" ] || [ "$TEST_MODE" = "directory" ]; then
    echo "  - test_output_directory/    Directory test results"
fi
echo ""
echo "Key validations performed:"
echo "  ✓ Single .txt file input"
echo "  ✓ Directory with multiple .txt files"
echo "  ✓ File merging and deduplication"
echo "  ✓ Filtering parameters passed to reports"
echo "  ✓ Comment and blank line filtering"
echo ""
echo "Usage examples:"
echo "  bash tests/run_test.sh                  # Run all tests"
echo "  bash tests/run_test.sh --single         # Run single file test only"
echo "  bash tests/run_test.sh --directory      # Run directory test only"
echo "  bash tests/run_test.sh --with-thermo    # Run with thermodynamic analysis"
echo "  bash tests/run_test.sh --all --with-thermo  # All tests with thermo"
echo ""
