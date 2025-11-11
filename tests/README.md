# Test Suite for Genome Primer Search Pipeline

This directory contains comprehensive tests for the Genome Primer Search pipeline, validating both core functionality and new features.

## Test Scripts

### 1. `run_test.sh` - Main Test Runner
Comprehensive test runner that orchestrates all test scenarios.

**Usage:**
```bash
# Run all tests
bash tests/run_test.sh

# Run specific test
bash tests/run_test.sh --single       # Single file test only
bash tests/run_test.sh --directory    # Directory test only

# Enable thermodynamic analysis
bash tests/run_test.sh --with-thermo
bash tests/run_test.sh --all --with-thermo
```

### 2. `test_single_file.sh` - Single File Input Test
Tests the pipeline with a single `.txt` file containing genome names/accessions.

**What it tests:**
- Single file input processing
- Parameter passing (min/max amplicon size)
- Deduplication within a single file
- Report generation with correct filter criteria

**Output:** `test_output_single/`

### 3. `test_directory.sh` - Directory Input Test
Tests the pipeline with a directory containing multiple `.txt` files.

**What it tests:**
- Directory input detection
- Multiple file merging
- Cross-file deduplication
- Comment and blank line filtering
- Parameter passing to report generation

**Output:** `test_output_directory/`

**Expected behavior:**
- Merges 3 files containing 6 total genome entries
- Deduplicates to 4 unique genomes
- Processes exactly 4 genomes

## Test Data

### `data/test_genomes.txt`
Single file with 2 test genomes:
- `GCA_000005845.2` (Mycoplasma genitalium, ~580kb)
- `Escherichia coli` (K-12, ~4.6Mb)

### `data/test_genomes_dir/`
Directory with 3 genome list files for testing merging and deduplication:

- **`genomes_bacteria.txt`**: 2 genomes (Mycoplasma, E. coli)
- **`genomes_additional.txt`**: 2 genomes + 2 duplicates (E. coli, GCA_000005845.2, Haemophilus influenzae)
- **`genomes_pathogens.txt`**: 1 genome (Streptococcus pneumoniae) + comments and blank lines

**Total entries:** 6
**Unique genomes:** 4 (after deduplication)

### `data/test_primers.txt`
Universal bacterial 16S rRNA primers:
- `16S_V3V4` (Illumina V3-V4 region)
- `16S_27F` (Universal forward primer)

### `data/test_references.fasta`
Reference sequences for thermodynamic analysis validation.

## Key Features Tested

### ✅ Input Handling
- [x] Single `.txt` file input
- [x] Directory containing multiple `.txt` files
- [x] File merging from directory
- [x] Deduplication (within and across files)
- [x] Comment line filtering (lines starting with `#`)
- [x] Blank line filtering

### ✅ Parameter Passing
- [x] `--min-size` parameter passed to report scripts
- [x] `--max-size` parameter passed to report scripts
- [x] Filter criteria displayed correctly in HTML reports
- [x] Filter criteria displayed correctly in TSV reports

### ✅ Pipeline Functionality
- [x] Genome downloading (from NCBI)
- [x] Primersearch execution
- [x] Amplicon filtering
- [x] Report generation
- [x] Thermodynamic analysis (optional)

### ✅ Output Validation
- [x] Correct number of genomes processed
- [x] No duplicate genome processing
- [x] HTML reports show actual filter parameters
- [x] TSV reports include filter criteria column

## Running Tests

### Quick Test (All Tests)
```bash
cd /path/to/genome-primer-search
bash tests/run_test.sh
```

### Individual Tests
```bash
# Test single file input
bash tests/test_single_file.sh

# Test directory input with deduplication
bash tests/test_directory.sh
```

### With Thermodynamic Analysis
```bash
# All tests with thermo
bash tests/run_test.sh --with-thermo

# Single test with thermo
bash tests/test_single_file.sh --with-thermo
```

## Expected Results

### Single File Test
- **Genomes processed:** 2
- **Filter criteria in reports:** "50 - 2,000 bp"
- **Output directory:** `test_output_single/`

### Directory Test
- **Files merged:** 3 (from `test_genomes_dir/`)
- **Total entries:** 6
- **Unique genomes processed:** 4 (after deduplication)
- **Filter criteria in reports:** "50 - 2,000 bp"
- **Output directory:** `test_output_directory/`

**Deduplication verification:**
- `Escherichia coli` appears in 2 files → processed once
- `GCA_000005845.2` appears in 2 files → processed once
- `GCA_000027305.1` appears once → processed once
- `GCA_000006885.1` appears once → processed once

## Troubleshooting

### Test fails with "file not found"
Ensure you're running from the repository root:
```bash
cd /path/to/genome-primer-search
bash tests/run_test.sh
```

### Genomes fail to download
- Check internet connection
- Verify NCBI datasets/efetch tools are installed
- Check NCBI servers are accessible

### Wrong number of genomes in directory test
- Verify all 3 `.txt` files exist in `tests/data/test_genomes_dir/`
- Check pipeline logs for deduplication messages
- Count genome files in `test_output_directory/results/genomes/`

## Adding New Tests

To add a new test scenario:

1. Create test data in `tests/data/`
2. Create new test script `tests/test_<name>.sh`
3. Add test to `run_test.sh` main runner
4. Update this README

## Test Coverage

| Feature | Single File Test | Directory Test |
|---------|-----------------|----------------|
| Single file input | ✅ | - |
| Directory input | - | ✅ |
| Deduplication | ✅ | ✅ |
| Parameter passing | ✅ | ✅ |
| Report generation | ✅ | ✅ |
| File merging | - | ✅ |
| Comment filtering | ✅ | ✅ |
| Blank line filtering | ✅ | ✅ |
| Thermo analysis | ✅ | ✅ |

## Notes

- Tests use small genomes to minimize download time and computational resources
- All tests can be run with `-resume` for faster iteration during development
- Test outputs are isolated to separate directories to avoid conflicts
- Color-coded output helps identify test results quickly
