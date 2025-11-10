# Genome Primer Search Pipeline

A Nextflow pipeline for downloading sequences from NCBI and searching for primer binding sites using EMBOSS primersearch.

## Quick Start

```bash
# Basic usage
nextflow run main.nf

# With custom inputs
nextflow run main.nf --genomes my_genomes.txt --primers my_primers.txt

# Resume after interruption
nextflow run main.nf -resume
```

## Features

- **Multiple input types**: Taxon names, genome assemblies, or individual nucleotide sequences
- **Smart filtering**: Removes circular genome artifacts (>10kb amplicons)
- **Scalable reports**: Long-format TSV works with any number of primers
- **Parallel processing**: Configurable concurrent downloads and searches
- **Modern architecture**: Nextflow DSL2 with modular processes

## Prerequisites

- **Nextflow** (>=22.10.0): `curl -s https://get.nextflow.io | bash`
- **Conda**: For automatic dependency management

The pipeline automatically installs: NCBI datasets CLI, NCBI E-utilities, and EMBOSS.

## Input Formats

### Genome/Sequence List

Create a text file with one entry per line. Mix any of these formats:

```
# Organism names
Escherichia coli
Homo sapiens

# Genome assembly accessions
GCA_000005845.2
GCF_000001405.40

# Nucleotide sequence accessions
U00096.3
NC_000913.3
MN986463.1

# Comments and empty lines OK
# Lines starting with # are ignored
```

**Supported accession formats:**
- **Assembly accessions**: `GCA_*` or `GCF_*` (downloads full genome assembly)
- **INSDC nucleotide**: 1-2 letters + 5-8 digits (e.g., U00096, MN986463, AF123456)
- **RefSeq nucleotide**: 2 letters + underscore (e.g., NC_000913, NM_001234)
- Version numbers supported: `.1`, `.2`, etc.

### Primers File

EMBOSS primersearch format (tab or space separated):

```
# Primer_Name Forward_Sequence Reverse_Sequence
16S_27F AGAGTTTGATCMTGGCTCAG TACGGYTACCTTGTTACGACTT
16S_V3V4 CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC
```

Supports IUPAC ambiguity codes (R, Y, W, M, K, etc.).

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genomes` | `genomes/genomes.txt` | Path to genome/sequence list |
| `--primers` | `primers/primers.txt` | Path to primers file |
| `--outdir` | `_output` | Output directory |
| `--mismatch` | `0` | Allowed mismatch percentage (0-100) |
| `--max_amplicon_size` | `10000` | Maximum amplicon size (bp) |
| `--min_amplicon_size` | `50` | Minimum amplicon size (bp) |
| `--max_downloads` | `4` | Concurrent genome downloads |
| `--max_searches` | `8` | Concurrent primersearch jobs |

## Output Structure

```
_output/
├── results/
│   ├── genomes/          # Downloaded FASTA files
│   ├── primersearch/     # Raw primersearch output
│   ├── filtered/         # Filtered results
│   └── filter_stats/     # JSON statistics per genome
└── reports/
    ├── summary.tsv       # Long-format summary (scales to any # of primers)
    ├── summary.html      # Interactive HTML report
    ├── amplicon_stats.tsv
    └── primer_stats.tsv
```

### Report Formats

**summary.tsv** - Long/normalized format for scalability:
```
Genome	Primer	Total_Found	Kept	Filtered_Circular	Filtered_Small	Mean_Size	Size_Range
E_coli	_GENOME_TOTAL_	36	14	22	0	-	-
E_coli	16S_27F	20	8	12	0	1375	465-1506
E_coli	16S_V3V4	15	6	9	0	465	465-465
```

Rows with `_GENOME_TOTAL_` provide per-genome aggregates. Easy to filter/analyze in Excel, R, or Python.

**primer_stats.tsv** - Per-primer summary:
```
Primer	Kept	Filtered_Total	Filtered_Circular	Filtered_Small	Genomes_Hit	Mean_Size	Size_Range
16S_27F	16	22	22	0	2	1375	465-1506
```

## Pipeline Workflow

```
Input → Download (datasets/efetch) → Extract → Primersearch → Filter → Reports
```

1. **Download**: Automatically detects input type and uses appropriate tool
2. **Extract**: Handles both zip archives and FASTA files
3. **Primersearch**: EMBOSS primersearch with configurable mismatches
4. **Filter**: Removes circular artifacts and size outliers
5. **Reports**: Generates scalable TSV and HTML reports

## Examples

### Test with small genomes
```bash
bash tests/run_test.sh
```

### Screen 16S primers against gut bacteria
```bash
nextflow run main.nf \
    --genomes genomes/genomes_CRC.txt \
    --primers primers/primers_16S.txt \
    --mismatch 5
```

### Download specific sequences
```bash
echo -e "NC_000913.3\nU00096.3" > my_sequences.txt
nextflow run main.nf --genomes my_sequences.txt
```

### High-throughput mode
```bash
nextflow run main.nf \
    --max_downloads 8 \
    --max_searches 16
```

## Filtering

The pipeline automatically filters amplicons:

- **Circular genome artifacts**: Amplicons >10kb (configurable with `--max_amplicon_size`)
- **Too small**: Amplicons <50bp (configurable with `--min_amplicon_size`)

Bacterial genomes are circular, so primers binding in opposite orientations can appear to span the entire genome. These false positives are automatically removed.

**Example**: E. coli with 16S primers
- Raw results: 36 amplicons found
- After filtering: 14 retained (22 circular artifacts removed)

Check `results/filter_stats/*.json` for detailed filtering information per genome.

## Troubleshooting

**No genome file found:**
- Check organism name spelling
- Try using assembly/nucleotide accession instead
- Verify organism has sequences in NCBI
- Check internet connectivity

**Pipeline hangs:**
- Check NCBI API status
- Reduce `--max_downloads` to avoid overwhelming NCBI
- Use `-resume` to restart from last successful step

**Conda environments slow to create:**
- First run creates environments (slow)
- Subsequent runs reuse environments (fast)
- Environments cached in `_work/conda/`

**Large amplicons (>100kb):**
- These are circular genome artifacts (now automatically filtered)
- Adjust threshold with `--max_amplicon_size`

## Architecture

Modern Nextflow DSL2 with modular processes:

```
main.nf                   # Entry point
workflows/primersearch.nf # Main workflow logic
modules/local/            # Individual process modules
  ├── download_genome.nf
  ├── extract_genome.nf
  ├── run_primersearch.nf
  ├── filter_amplicons.nf
  └── generate_reports.nf
bin/                      # Helper scripts
  ├── filter_amplicons.py
  └── generate_html_report.py
```

Easy to extend with new processes or modify existing ones.

## Performance Tips

1. **Parallel execution**: Adjust `--max_downloads` and `--max_searches` based on your system
2. **Resume capability**: Always use `-resume` for interrupted runs
3. **Test first**: Use `tests/run_test.sh` with small genomes before large datasets
4. **Work directory**: Clean old runs with `rm -rf _work/` to free space
5. **Scalability**: Pipeline handles 10-1000+ genomes with proper settings

## Citation

If you use this pipeline, please cite:
- EMBOSS: Rice et al. (2000) EMBOSS: The European Molecular Biology Open Software Suite
- NCBI Datasets: https://www.ncbi.nlm.nih.gov/datasets/

## License

MIT License
