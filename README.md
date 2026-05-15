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
- **Thermodynamic analysis**: Optional off-target amplification risk assessment using nearest-neighbor Tm calculations
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

### Basic Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genomes` | `data/genomes/` | Path to a genome/sequence list `.txt` file or a directory of `.txt` files |
| `--primers` | _required_ | Path to primers file (e.g. `data/primers/16S_primers.txt`) |
| `--outdir` | `_output` | Output directory |
| `--mismatch` | `0` | Allowed mismatch percentage (0-100) |
| `--max_primer_size` | `100` | Maximum primer length (bases) |
| `--max_amplicon_size` | `10000` | Maximum amplicon size (bp) |
| `--min_amplicon_size` | `50` | Minimum amplicon size (bp) |
| `--max_downloads` | `4` | Concurrent genome downloads |
| `--max_searches` | `8` | Concurrent primersearch jobs |

### Thermodynamic Analysis Parameters (Optional)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--enable_thermo_analysis` | `false` | Enable thermodynamic off-target analysis |
| `--thermo_references` | `references/target_sequences.fasta` | Reference target sequences (FASTA) |
| `--thermo_master_mix` | `DreamTaq` | PCR master mix preset (see pcr_conditions/master_mixes.json) |
| `--thermo_annealing_temp` | `60` | Annealing temperature (°C) |
| `--thermo_delta_tm_high` | `5` | ΔTm threshold for high off-target risk (°C) |
| `--thermo_delta_tm_medium` | `10` | ΔTm threshold for medium off-target risk (°C) |

**Available master mixes:** DreamTaq, Q5_HighFidelity, Phusion_HF, Phusion_GC, OneTaq, GoTaq, KAPA_HiFi, Platinum_Taq, Custom

## Output Structure

```
_output/
├── results/
│   ├── genomes/          # Downloaded FASTA files
│   ├── primersearch/     # Raw primersearch output
│   ├── filtered/         # Filtered results
│   └── filter_stats/     # JSON statistics per genome
├── reports/
│   ├── summary.tsv       # Long-format summary (scales to any # of primers)
│   ├── summary.html      # Interactive HTML report
│   ├── amplicon_stats.tsv
│   ├── primer_stats.tsv
│   └── failed_genomes.tsv  # Genomes that failed to download — written only if any failed
└── thermo_reports/       # Optional: Thermodynamic analysis (if enabled)
    ├── thermo_analysis.tsv       # Per-amplicon thermodynamic metrics
    ├── offtarget_summary.tsv     # High-risk off-target amplifications
    ├── primer_specificity.tsv    # Per-primer specificity scores
    └── thermo_analysis.html      # Interactive thermodynamic report
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
    --genomes data/genomes/human_gut_microbiome_CRC_associated.txt \
    --primers data/primers/16S_primers.txt \
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

## Thermodynamic Off-Target Analysis

The pipeline includes an optional thermodynamic analysis module that estimates the probability of off-target amplification using melting temperature (Tm) calculations.

### When to Use

Enable thermodynamic analysis when:
- Screening primers against genomes with potential off-target sequences (e.g., host contamination)
- Designing multiplexed assays where specificity is critical
- Evaluating primer performance across multiple organisms
- Assessing risk of cross-reactivity with related species

### How It Works

For each primer-target pair, the pipeline calculates:

1. **Target Tm**: Melting temperature of the primer binding to the intended reference target
2. **Actual Tm**: Melting temperature accounting for mismatches in the observed binding
3. **ΔTm**: Difference between target and actual Tm (ΔTm = Tm_target - Tm_actual)
4. **Mismatch profile**: Number and estimated location of mismatches (3′ end critical)
5. **Amplification probability**: Risk classification (high/medium/low) based on decision rules:
   - **High risk**: ΔTm < 5°C and minimal 3′ mismatches
   - **Medium risk**: ΔTm 5-10°C or some 3′ mismatches
   - **Low risk**: ΔTm > 10°C or significant 3′ mismatches

### Classification Logic: Intended Targets vs Off-Targets

The analysis distinguishes between **intended targets** (expected amplification) and **off-targets** (unintended amplification):

- **Intended Target**: An amplicon with zero mismatches in both primers that matches a reference sequence
  - Classification: `target` (not a risk)
  - P(Amplification): 1.0 (100% expected)
  - ΔTm: 0°C (perfect match)
  - These are your desired amplicons - the primers are working as designed

- **Off-Target**: Any amplicon that is NOT an intended target
  - Classification: `high`, `medium`, or `low` risk
  - P(Amplification): 0.10 - 0.85 depending on ΔTm and mismatches
  - ΔTm: > 0°C (indicates deviation from intended target)
  - These amplicons may compete with or interfere with your intended targets

**Important**: A perfect match (0 mismatches) to a genome sequence is **only** classified as an intended target if a matching reference sequence is provided. Without a reference, all amplicons are treated as potential off-targets.

**Example**:
- Primer `16S_V3V4` amplifies E. coli 16S gene with 0 mismatches
- Reference `16S_V3V4_reference` exists → Classified as `target` (good!)
- No reference provided → Classified as `high` risk (needs review)

This distinction ensures specificity metrics accurately reflect primer performance: high-specificity primers amplify their intended targets without generating high-risk off-targets.

### Tm Calculations

Uses BioPython's nearest-neighbor thermodynamics with salt corrections (Owczarzy et al. 2008):
- Accounts for buffer composition (Na⁺, Mg²⁺, dNTPs)
- Applies penalties for mismatches (~1.5°C per mismatch)
- Supports PCR master mix presets with validated buffer compositions

### Usage Example

```bash
# Enable thermodynamic analysis with DreamTaq buffer
nextflow run main.nf \
    --genomes data/genomes/human_gut_microbiome_common.txt \
    --primers data/primers/16S_primers.txt \
    --enable_thermo_analysis \
    --thermo_references data/references/16S_ref.fasta \
    --thermo_master_mix DreamTaq \
    --thermo_annealing_temp 60

# Use Q5 high-fidelity polymerase with higher annealing temp
nextflow run main.nf \
    --enable_thermo_analysis \
    --thermo_master_mix Q5_HighFidelity \
    --thermo_annealing_temp 65

# Adjust sensitivity thresholds
nextflow run main.nf \
    --enable_thermo_analysis \
    --thermo_delta_tm_high 3 \
    --thermo_delta_tm_medium 8
```

### Reference Sequences

Reference sequences define the **intended target** for each primer pair. They're used to calculate ΔTm between perfect binding and actual binding.

**How Reference Matching Works:**

1. **Primer name** (from primer file): `fadA`
2. **Reference ID** (from FASTA file): `fadA_reference`
3. **Matching**: Primer name must match everything before `_reference` in the FASTA ID

**You can provide:**
- **All references** (one FASTA file with all target genes)
- **Some references** (only the genes you care about)
- **No references** (all primers analyzed as off-targets)

Missing references generate warnings but analysis continues.

**Example 1: Multiple Virulence Genes**

Your primer file (`primers/virulence_genes.txt`):
```
fadA ATGAAACGCATCGCACAGCT CGTAGCAGTAGCATTCGTAG
fap2 GCTAGCTTAGGCCTAGCTAG AATTGGCCAATTGGCCAATT
radA GGCCTAGGCCTAGGCCTAGG CCAATTCCAATTCCAATTCC
```

Your reference file (`references/virulence_targets.fasta`):
```
>fadA_reference F. nucleatum fadA adhesin gene
ATGAAACGCATCGCACAGCT...full gene sequence...CGTAGCAGTAGCATTCGTAG
>fap2_reference F. nucleatum fap2 adhesin gene
GCTAGCTTAGGCCTAGCTAG...full gene sequence...AATTGGCCAATTGGCCAATT
>radA_reference DNA repair radA gene
GGCCTAGGCCTAGGCCTAGG...full gene sequence...CCAATTCCAATTCCAATTCC
```

Each primer is analyzed independently. You can include only the genes you need to validate.

**Example 2: 16S Primers**

```
>16S_V3V4_Novogene_reference E.coli 16S V3-V4 region
CCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCT...
>16S_FullLength_reference E.coli 16S full length
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACAC...
```

**Important:** If you're screening primers against multiple genomes but only care about specific target genes, provide ONLY those genes as references. Off-target hits in other genomes will be flagged appropriately.

### PCR Master Mix Configuration

Master mix presets are defined in `pcr_conditions/master_mixes.json`. Each preset includes:
- **Buffer composition** (Na⁺, Mg²⁺, dNTPs, Tris, KCl) - Used for Tm calculations
- Polymerase type and fidelity
- Typical annealing/extension temperatures (for reference)

**Common presets:**
- **DreamTaq**: Standard Taq, 2mM Mg²⁺, suitable for routine PCR
- **Q5_HighFidelity**: High-fidelity polymerase, higher annealing temps
- **Phusion_HF/GC**: Ultra high-fidelity, optimized for GC-rich templates
- **Custom**: Template for user-defined conditions

**About `--thermo_annealing_temp`:**

This parameter records the **actual annealing temperature you use in your PCR**, not necessarily the optimal temperature. It's stored in the report metadata but doesn't affect Tm calculations.

- **Tm calculations** use buffer composition (from master mix preset)
- **Annealing temp** is what YOU set on your thermocycler
- They can differ: You might use 60°C annealing even if calculated Tm is 65°C

**Example:**
```bash
# DreamTaq has typical annealing temp of 60°C, but you use 58°C
nextflow run main.nf \
    --enable_thermo_analysis \
    --thermo_master_mix DreamTaq \
    --thermo_annealing_temp 58  # Your actual PCR condition
```

The analysis will use DreamTaq's buffer chemistry for Tm calculation, but record that you used 58°C for annealing.

### Output Reports

**thermo_analysis.tsv** - Detailed per-amplicon metrics:
```
Primer	Genome	Amplicon_Length	Risk_Classification	P_Amplification	Delta_Tm_Average	Total_Mismatches
16S_V3V4	E_coli	465	low	0.010	15.2	0
16S_V3V4	H_sapiens_MT	470	high	0.723	3.1	1
```

**offtarget_summary.tsv** - Flagged high/medium risk off-targets:
```
Primer	Genome	Risk_Classification	Delta_Tm_Average	Recommendation
16S_V3V4	H_sapiens_MT	HIGH	3.1	REVIEW REQUIRED - High amplification probability
```

**primer_specificity.tsv** - Per-primer performance:
```
Primer	Total_Hits	High_Risk_OffTargets	Specificity_Score	Specificity_Rating
16S_V3V4	25	2	92.0	Excellent
16S_V4_EMP	18	5	72.2	Good
```

### Interpreting Results

- **Specificity Score**: 100 = perfect (no high-risk off-targets), 0 = poor
- **High-risk off-targets**: Should be reviewed; may require primer redesign
- **ΔTm < 5°C**: Off-target likely to amplify under standard conditions
- **Perfect match (0 mismatches) to off-target**: Strong cross-reactivity risk

### Understanding Alignment Results

The thermodynamic report includes detailed primer-genome alignments with the following features:

#### Degenerate Primers and Variant Selection

**What are degenerate primers?** Primers containing IUPAC ambiguity codes like `[TCY]` (T or C or Y) represent multiple possible sequences. For example, `AG[TCY]GTT` represents 3 different primers: `AGTGTT`, `AGCGTT`, and `AGYGTT`.

**How alignment works:**
1. EMBOSS primersearch matches primers using **all possible variants** (IUPAC-aware matching)
2. Our alignment tool selects the **best matching variant** for visualization
3. For each degenerate position `[ABC]`, the base that matches the genome is selected

**Example:**
```
Degenerate primer: AG[CAM]GTT[TCY]GAT
Genome sequence:   AGAGTTTGAT
Matched variant:   AGAGTTTGAT  (selects A from [CAM], T from [TCY])
```

#### Two Mismatch Counts

You'll see two different mismatch counts displayed:

1. **Variant mismatches** (shown in alignment): Mismatches between the selected variant and the genome
2. **Degenerate mismatches** (primersearch count): Mismatches considering all possible variants

**When they differ:**
```
Forward: 4 variant mismatches | 0 degenerate mismatches (1 at 3' end)
```

This means:
- The selected variant has 4 mismatches with the genome
- But another variant exists (from the degenerate code) that matches perfectly
- Primersearch found 0 mismatches because it tries all variants

**Why this matters:** For **intended targets** with 0 degenerate mismatches, this indicates the primer is working correctly - one of its variants matches perfectly. For **off-targets**, having many variant mismatches but few degenerate mismatches suggests the degenerate codes are helping the primer bind broadly.

#### Interpreting Intended Target Results

**Expected for perfect intended targets:**
- Degenerate mismatches: 0
- Variant mismatches: 0 (if all degenerate positions happen to match)
- ΔTm: 0°C
- P(Amplification): 1.0

**If you see mismatches on intended targets:**
- Check the "Matched Variant" section to see which specific variant was used
- Non-zero variant mismatches with 0 degenerate mismatches is expected when primers contain ambiguity codes
- The alignment shows one specific variant path through the degenerate code space

#### 3' End Mismatches

Mismatches within the **last 5 bases** from the 3' end are critical for PCR:
- DNA polymerase extends from the 3' end
- Even 1 mismatch at the 3' end significantly reduces amplification efficiency
- Mismatches at position 0-1 from 3' end: ~90% reduction in efficiency
- Mismatches at position 2-4 from 3' end: ~50% reduction in efficiency
- Mismatches at 5+ from 3' end: Minimal effect

The report highlights these as they're more important than 5' mismatches.

#### Gaps in Alignment

**What gaps mean:**
```
CCTACG--GGCGGCT-GCAG
.|  ||  ||| .|| ||||
GC--CGGTGGC-ACTGGCAG
Warning: 6 gaps in alignment
```

Multiple gaps (>3) indicate:
- Position extraction may be imprecise
- Could occur at sequence boundaries
- Complex primers with repetitive regions
- The genomic context doesn't allow perfect alignment

**This doesn't mean the result is wrong** - primersearch found the binding site. But the alignment is showing you that visualizing it as a linear match is challenging.

#### Color Coding

In the HTML report:
- **Green**: 0 mismatches (perfect match)
- **Yellow**: 1-2 mismatches (acceptable)
- **Red**: 3+ mismatches (potentially problematic)

These colors apply to the variant mismatch count shown in the alignment.

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
main.nf                           # Entry point
workflows/
  ├── primersearch.nf             # Main workflow
  └── thermo_analysis.nf          # Optional thermodynamic analysis
modules/local/                    # Process modules
  ├── download_genome.nf
  ├── extract_genome.nf
  ├── run_primersearch.nf
  ├── filter_amplicons.nf
  ├── generate_reports.nf
  ├── thermodynamic_analysis.nf   # Tm calculations
  └── generate_thermo_reports.nf  # Thermodynamic reports
bin/                              # Helper scripts
  ├── filter_amplicons.py
  ├── generate_html_report.py
  ├── thermodynamic_analysis.py   # Tm calculation engine
  └── generate_thermo_report.py   # Report generation
pcr_conditions/
  └── master_mixes.json           # PCR buffer presets
references/
  └── target_sequences.fasta      # Reference target sequences
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
- BioPython: Cock et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics
- Thermodynamics: Owczarzy et al. (2008) Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations

## License

MIT License
