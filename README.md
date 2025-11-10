# Genome Primer Search Pipeline

A Nextflow pipeline for downloading reference genomes from NCBI and searching for primer binding sites using EMBOSS primersearch.

## Overview

This pipeline automates the process of:
1. Downloading reference genomes from NCBI (supports both organism names and assembly accessions)
2. Extracting genome sequences
3. Running primersearch to find primer binding sites
4. Generating a summary report of all results

## Features

- **Flexible Input**: Supports both organism taxon names (e.g., "Escherichia coli") and assembly accession IDs (e.g., "GCA_000005845.2")
- **Universal Scope**: Works with any genome type (bacterial, archaeal, fungal, human, etc.)
- **16S rRNA Support**: Includes universal bacterial 16S primers for microbiome studies
- **Parallel Processing**: Downloads and analyzes multiple genomes simultaneously
- **Automatic Retry**: Retries failed downloads up to 2 times

## Prerequisites

### Required Software

- **Nextflow** (>=22.10.0)
  ```bash
  curl -s https://get.nextflow.io | bash
  ```

- **Conda** (recommended for automatic dependency management)
  ```bash
  # Install Miniconda
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ```

The pipeline will automatically install:
- **NCBI datasets CLI** (via conda-forge)
- **EMBOSS** (for primersearch, via bioconda)

## Input Files

### 1. Genome List (`genomes.txt`)

A text file supporting two input formats:

**Format 1: Organism taxon names** (scientific names)
```
Homo sapiens
Escherichia coli
Fusobacterium nucleatum
```

**Format 2: Assembly accession IDs** (GCA_* or GCF_*)
```
GCA_000005845.2
GCF_000001405.40
```

**Mixed format** (both types in same file)
```
# Human genome
Homo sapiens

# CRC-associated bacteria
Fusobacterium nucleatum
Bacteroides fragilis

# Small test genome (by accession)
GCA_000005845.2
```

**Features**:
- Lines starting with `#` are treated as comments
- Empty lines are ignored
- Mix taxon names and accessions freely
- Bracketed names supported (e.g., `[Clostridium] symbiosum`)

The provided `genomes.txt` contains:
- Human genome
- 17 CRC-associated bacterial species
- 1 test genome (Mycoplasma genitalium, ~580kb)

### 2. Primers File (`primers.txt` or `primers_16S.txt`)

EMBOSS primersearch format:

```
# Primer_Name Forward_Sequence Reverse_Sequence
PrimerA ATGCGATCGATCG CGTAGCTAGCTAG
PrimerB GCTAGCTAGCTAG ATCGATCGATCGA
```

**Provided primer files**:
- **`primers.txt`**: Gene-specific primers (FadA, Fap2, RadD, Human β-globulin)
- **`primers_16S.txt`**: Universal bacterial 16S rRNA primers
  - 27F/1492R: Full-length 16S (~1,400 bp)
  - 341F/785R: V3-V4 regions (~450 bp)
  - 515F/806R: V4 region (~291 bp, Earth Microbiome Project standard)

## Usage

### Basic Usage

Run with default parameters (uses `genomes.txt` and `primers.txt`):

```bash
nextflow run main.nf
```

### Custom Input Files

Specify custom genome list and/or primers file:

```bash
# Use 16S primers instead of gene-specific primers
nextflow run main.nf --primers primers_16S.txt

# Use custom genome list
nextflow run main.nf --genomes my_genomes.txt

# Use both custom files
nextflow run main.nf --genomes my_genomes.txt --primers my_primers.txt
```

### Additional Parameters

```bash
nextflow run main.nf \
    --genomes my_genomes.txt \
    --primers primers_16S.txt \
    --outdir my_results \
    --mismatch 10
```

### Using Conda (Recommended)

The pipeline automatically uses conda for dependencies:

```bash
nextflow run main.nf --primers primers_16S.txt
```

### Resume Failed Run

```bash
nextflow run main.nf -resume
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genomes` | `genomes.txt` | Path to file with genome names/accessions |
| `--primers` | `primers.txt` | Path to primers file (EMBOSS format) |
| `--outdir` | Auto-generated | Output directory (auto: `_output_TIMESTAMP_PRIMERNAME`) |
| `--mismatch` | `0` | Allowed mismatch percentage (0-100) |

**Note on output directory**: By default, the pipeline creates a timestamped output directory based on the run time and primer filename used (e.g., `_output_2025-11-10_20-45-30_primers`). This allows you to run the pipeline multiple times with different primers without overwriting results. You can override this with `--outdir` if needed.

## Input Format Details

### Organism Names vs. Accession IDs

The pipeline automatically detects input type:

**Organism Names**:
- Downloads **reference genome** for the taxon
- Example: `Escherichia coli` → downloads reference E. coli genome
- Use for: Getting representative genomes

**Accession IDs**:
- Downloads **specific assembly** by accession
- Example: `GCA_000005845.2` → downloads exact assembly
- Use for: Reproducible research, specific genome versions

**When to use accessions**:
- Organism has multiple strains/assemblies
- Need specific genome version
- Taxonomic names are ambiguous or reclassified (e.g., `[Clostridium] symbiosum`)

### Special Cases

**Bracketed taxonomy names**: Some organisms have uncertain taxonomy and use brackets:
```
[Clostridium] symbiosum
[Eubacterium] rectale
```

**Human genome**: Use either:
```
Homo sapiens          # Gets latest reference (GRCh38)
GCF_000001405.40      # Specific version (GRCh38.p14)
```

## Output Structure

Each run creates a timestamped output directory named with the format:
`_output_YYYY-MM-DD_HH-MM-SS_PRIMERNAME`

Example directory structure:

```
_output_2025-11-10_20-45-30_primers/
├── genome/
│   ├── Homo_sapiens_genome.fasta
│   ├── Escherichia_coli_genome.fasta
│   ├── GCA_000005845.2_genome.fasta
│   └── ...
├── primersearch_results/
│   ├── Homo_sapiens_primersearch.out
│   ├── Escherichia_coli_primersearch.out
│   ├── GCA_000005845.2_primersearch.out
│   └── ...
├── summary_report.txt
└── reports/
    ├── report.html
    ├── timeline.html
    ├── trace.txt
    └── flowchart.html

_output_2025-11-10_21-30-15_primers_16S/
├── genome/
├── primersearch_results/
├── summary_report.txt
└── reports/
```

**Benefits**:
- Each run with different primers creates a separate output directory
- Previous results are never overwritten
- Easy to compare results from different primer sets
- Directory name shows when the run was performed and which primers were used

### Output Files

- **genome/**: Downloaded and extracted genome FASTA files
- **primersearch_results/**: Individual primersearch results for each genome
- **summary_report.txt**: Aggregated summary of all primer hits across genomes
- **reports/report.html**: Detailed execution report with resource usage
- **reports/timeline.html**: Visual timeline of pipeline execution
- **reports/trace.txt**: Process execution trace for debugging
- **reports/flowchart.html**: Pipeline DAG visualization

## Pipeline Workflow

```
┌─────────────┐
│ genomes.txt │
└──────┬──────┘
       │
       ▼
┌──────────────────┐
│ DOWNLOAD_GENOME  │ ◄── Detects taxon vs accession
│                  │     Parallel processing
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ EXTRACT_GENOME   │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐     ┌──────────────┐
│ RUN_PRIMERSEARCH │ ◄───│ primers.txt  │
└────────┬─────────┘     └──────────────┘
         │
         ▼
┌──────────────────┐
│ COLLECT_RESULTS  │
└────────┬─────────┘
         │
         ▼
    ┌─────────┐
    │ Results │
    └─────────┘
```

## Examples

### Example 1: Screen 16S primers against gut bacteria

```bash
nextflow run main.nf \
    --genomes genomes.txt \
    --primers primers_16S.txt \
    --mismatch 5
```

### Example 2: Add new genomes by accession

1. Edit `genomes.txt`:
   ```bash
   echo "GCA_000146045.2" >> genomes.txt  # Salmonella Typhimurium
   ```

2. Run pipeline:
   ```bash
   nextflow run main.nf -resume
   ```

### Example 3: Custom primer file for specific gene

1. Create custom primers in EMBOSS format:
   ```bash
   cat > my_primers.txt << 'EOF'
   # Gene_X_primers
   ATGCGATCGATCG CGTAGCTAGCTAG
   EOF
   ```

2. Run with custom primers:
   ```bash
   nextflow run main.nf --primers my_primers.txt
   ```

### Example 4: Human genome with custom primers

1. Create file with just human:
   ```bash
   echo "Homo sapiens" > human_only.txt
   ```

2. Run analysis:
   ```bash
   nextflow run main.nf \
       --genomes human_only.txt \
       --primers primers_16S.txt
   ```

## Troubleshooting

### Issue: "PackagesNotFoundError: ncbi-datasets-cli"

**Solution**: The conda channel was incorrect - this is now fixed (uses conda-forge instead of bioconda)

### Issue: "Clostridium symbiosum not found"

**Solution**: Use bracketed name `[Clostridium] symbiosum` or accession ID instead

### Issue: "No genome file found"

**Solution**:
- Check organism name spelling in `genomes.txt`
- Try using assembly accession instead of taxon name
- Verify the organism has a genome in NCBI
- Check internet connectivity
- Look at `.nextflow.log` for detailed error messages

### Issue: Pipeline hangs during download

**Solution**:
- Check NCBI API status
- Reduce number of genomes to download
- Use `-resume` to restart from last successful step
- Increase timeout in `nextflow.config` (process.withName:DOWNLOAD_GENOME.time)

### Issue: Conda environment creation slow

**Solution**:
- First run creates conda environments (slow)
- Subsequent runs reuse environments (fast)
- Environments cached in `work/conda/`

## Performance Tips

1. **Parallel Execution**: Pipeline automatically processes multiple genomes in parallel
2. **Resume Capability**: Use `-resume` to continue from last successful step
3. **Resource Allocation**: Adjust CPU/memory in `nextflow.config` based on your system
4. **Conda Caching**: Conda environments are cached after first creation
5. **Small Test Genome**: Use `GCA_000005845.2` (Mycoplasma, ~580kb) to test pipeline quickly

## 16S Primer Information

The included `primers_16S.txt` contains three primer sets for bacterial 16S rRNA amplification:

### 27F/1492R (Full-Length)
- **Target**: Nearly complete 16S gene (~1,400 bp)
- **Use case**: Species-level identification, phylogenetic analysis
- **Reference**: Weisburg et al. (1991)

### 341F/785R (V3-V4)
- **Target**: V3-V4 hypervariable regions (~450 bp)
- **Use case**: Microbiome diversity studies, highest OTU recovery
- **Best for**: Illumina MiSeq, most reproducible

### 515F/806R (V4)
- **Target**: V4 hypervariable region (~291 bp)
- **Use case**: Earth Microbiome Project standard
- **Note**: Also amplifies archaeal 16S

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319.
- **NCBI Datasets**: https://www.ncbi.nlm.nih.gov/datasets/
- **EMBOSS**: Rice, P., et al. (2000). EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics 16(6):276-7.

## License

MIT License

## Contact

For issues or questions, please open an issue on the repository.
