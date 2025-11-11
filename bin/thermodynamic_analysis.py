#!/usr/bin/env python3
"""
Thermodynamic Analysis of Primer Binding Sites

Analyzes primersearch results to calculate:
- Melting temperatures (Tm) for target and off-target amplicons
- Delta Tm between intended and actual targets
- Mismatch profiles
- Amplification probability estimates

Uses BioPython for Tm calculations with nearest-neighbor thermodynamics
"""

import argparse
import json
import re
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Align import PairwiseAligner
    import math
except ImportError:
    print("Error: BioPython is required. Install with: conda install -c conda-forge biopython", file=sys.stderr)
    sys.exit(1)

# Physical constants
R_CAL = 1.987  # Gas constant in cal/(mol·K)
R_KCAL = R_CAL / 1000.0  # Gas constant in kcal/(mol·K)


class ThermoAnalyzer:
    """Main class for thermodynamic analysis of primer binding"""

    def __init__(self, pcr_conditions: Dict, primers: Dict, references: Dict,
                 genome_fasta: Optional[Path] = None,
                 annealing_temp: float = 60.0, delta_tm_thresholds: Tuple[float, float] = (5.0, 10.0)):
        """
        Initialize thermodynamic analyzer

        Args:
            pcr_conditions: Dict with buffer composition (Na, Mg, dNTPs in mM)
            primers: Dict mapping primer names to (forward_seq, reverse_seq) tuples
            references: Dict mapping primer names to reference target sequences
            genome_fasta: Path to genome FASTA file for extracting binding sequences
            annealing_temp: Annealing temperature in Celsius
            delta_tm_thresholds: Tuple of (high_risk, medium_risk) ΔTm thresholds in °C
        """
        self.pcr_conditions = pcr_conditions
        self.primers = primers
        self.references = references
        self.genome_fasta = genome_fasta
        self.annealing_temp = annealing_temp
        self.delta_tm_high, self.delta_tm_medium = delta_tm_thresholds

    def calculate_tm(self, primer_seq: str, target_seq: str,
                     mismatch_count: int = 0) -> float:
        """
        Calculate melting temperature using nearest-neighbor method

        Args:
            primer_seq: Primer sequence
            target_seq: Target/template sequence (same orientation as primer)
            mismatch_count: Number of mismatches (affects calculation)

        Returns:
            Melting temperature in Celsius
        """
        try:
            # Use Wallace rule for very short primers
            if len(primer_seq) < 14:
                tm = self._calculate_tm_wallace(primer_seq)
            else:
                # Use nearest-neighbor method with salt corrections
                tm = mt.Tm_NN(
                    Seq(primer_seq),
                    Na=self.pcr_conditions.get('Na', 50.0),
                    Mg=self.pcr_conditions.get('Mg', 2.0),
                    dNTPs=self.pcr_conditions.get('dNTPs', 0.2),
                    saltcorr=7  # Owczarzy et al. 2008 salt correction
                )

            # Apply penalty for mismatches (empirical: ~1-2°C per mismatch)
            if mismatch_count > 0:
                tm -= (1.5 * mismatch_count)

            return round(tm, 2)

        except Exception as e:
            print(f"Warning: Tm calculation failed for {primer_seq[:20]}...: {e}", file=sys.stderr)
            return 0.0

    def _calculate_tm_wallace(self, seq: str) -> float:
        """Wallace rule for short primers: Tm = 2(A+T) + 4(G+C)"""
        at_count = seq.upper().count('A') + seq.upper().count('T')
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        return 2 * at_count + 4 * gc_count

    def calculate_binding_probability_from_tm(self, tm_actual: float, tm_target: float,
                                              temp_celsius: float) -> float:
        """
        Calculate binding probability using a simplified empirical model

        Based on the relationship between Tm, annealing temperature, and binding.
        At T_anneal = Tm, approximately 50% of primers are bound.

        Args:
            tm_actual: Actual Tm with mismatches (°C)
            tm_target: Perfect match Tm (°C)
            temp_celsius: Annealing temperature (°C)

        Returns:
            Binding probability (0-1)
        """
        # Use a sigmoid function based on the difference between Tm and annealing temp
        # When Tm >> T_anneal: high binding (P → 1)
        # When Tm ≈ T_anneal: moderate binding (P ≈ 0.5)
        # When Tm << T_anneal: low binding (P → 0)

        # The "steepness" parameter controls how sharply binding drops off
        # Typical value: 1 degree Tm change ≈ ~10% change in binding
        steepness = 0.15  # Empirical parameter

        # Calculate binding probability using sigmoid
        # P_bind = 1 / (1 + exp(-steepness * (Tm - T_anneal)))
        delta_t = tm_actual - temp_celsius

        try:
            p_bind = 1.0 / (1.0 + math.exp(-steepness * delta_t))
        except OverflowError:
            p_bind = 1.0 if delta_t > 0 else 0.0

        return p_bind

    def analyze_mismatch_profile(self, primer_seq: str, mismatches: int,
                                  primer_strand: str) -> Dict:
        """
        Analyze mismatch profile and assign risk score

        Args:
            primer_seq: Primer sequence
            mismatches: Number of mismatches reported by primersearch
            primer_strand: 'forward' or 'reverse'

        Returns:
            Dict with mismatch analysis
        """
        # Count 3' mismatches (last 3 bases are critical for extension)
        # Note: primersearch only reports total mismatches, not positions
        # This is a limitation - we assume uniform distribution

        primer_length = len(primer_seq)
        three_prime_region = 3  # Last 3 bases

        # Estimate probability that mismatches are in 3' region
        # Assumes random distribution (conservative estimate)
        prob_3prime_mismatch = min(1.0, (three_prime_region / primer_length) * mismatches)

        return {
            'total_mismatches': mismatches,
            'primer_length': primer_length,
            'estimated_3prime_mismatches': round(prob_3prime_mismatch, 2),
            '3prime_region_size': three_prime_region,
            'mismatch_density': round(mismatches / primer_length, 3) if primer_length > 0 else 0
        }

    def calculate_amplification_probability(self, primer_seq: str,
                                           tm_actual: float,
                                           tm_target: float,
                                           delta_tm: float,
                                           mismatch_profile: Dict,
                                           primer_strand: str,
                                           is_intended_target: bool = False) -> Tuple[float, str]:
        """
        Estimate probability of amplification using Tm-based empirical model

        Uses sigmoid function relating Tm and annealing temperature for binding,
        combined with extension efficiency based on 3' mismatches.

        Args:
            primer_seq: Primer sequence
            tm_actual: Actual Tm with mismatches (°C)
            tm_target: Perfect match target Tm (°C)
            delta_tm: ΔTm = Tm_target - Tm_actual
            mismatch_profile: Mismatch analysis dict
            primer_strand: 'forward' or 'reverse'
            is_intended_target: Whether this is the intended target (not an off-target)

        Returns:
            Tuple of (probability 0-1, classification 'target'/'high'/'medium'/'low')
        """
        # Intended targets should not be classified as "risk"
        # Perfect match to intended target = expected behavior
        if is_intended_target:
            return (1.0, 'target')

        # Get mismatch information
        estimated_3prime = mismatch_profile['estimated_3prime_mismatches']
        total_mm = mismatch_profile['total_mismatches']

        # Calculate binding probability based on Tm and annealing temperature
        # Uses sigmoid relationship
        p_bind = self.calculate_binding_probability_from_tm(
            tm_actual,
            tm_target,
            self.annealing_temp
        )

        # Calculate extension efficiency
        # 3' mismatches severely reduce polymerase extension
        # Use exponential penalty: P_extend = exp(-k * 3prime_mismatches)
        k_extension = 1.2  # Empirical constant for 3' mismatch penalty
        p_extend = math.exp(-k_extension * estimated_3prime)

        # Overall probability: must both bind AND extend
        p_amplification = p_bind * p_extend

        # Classify based on probability thresholds
        if p_amplification >= 0.7:
            classification = 'high'
        elif p_amplification >= 0.3:
            classification = 'medium'
        else:
            classification = 'low'

        return (p_amplification, classification)

    def analyze_amplicon(self, primer_name: str, genome_name: str,
                        forward_seq: str, reverse_seq: str,
                        forward_mm: int, reverse_mm: int,
                        amplicon_length: int,
                        forward_pos: Optional[int] = None,
                        reverse_pos: Optional[int] = None) -> Dict:
        """
        Complete thermodynamic analysis for one amplicon

        Args:
            primer_name: Name of primer pair
            genome_name: Genome/sequence name
            forward_seq: Forward primer sequence
            reverse_seq: Reverse primer sequence
            forward_mm: Forward primer mismatches
            reverse_mm: Reverse primer mismatches
            amplicon_length: Amplicon size in bp
            forward_pos: Forward primer genomic position (for alignment extraction)
            reverse_pos: Reverse primer genomic position (for alignment extraction)

        Returns:
            Dict with complete thermodynamic analysis
        """
        # Get reference sequence for this primer
        reference_seq = self.references.get(primer_name)

        if reference_seq is None:
            print(f"Warning: No reference sequence for primer {primer_name}", file=sys.stderr)
            is_target = False
            tm_target_fwd = tm_target_rev = 0.0
        else:
            # Check if this is the intended target (will refine this logic)
            is_target = (forward_mm == 0 and reverse_mm == 0)

            # Calculate Tm for perfect target binding
            tm_target_fwd = self.calculate_tm(forward_seq, reference_seq, mismatch_count=0)
            tm_target_rev = self.calculate_tm(reverse_seq, reference_seq, mismatch_count=0)

        # Calculate Tm for actual binding (with mismatches)
        tm_actual_fwd = self.calculate_tm(forward_seq, "", mismatch_count=forward_mm)
        tm_actual_rev = self.calculate_tm(reverse_seq, "", mismatch_count=reverse_mm)

        # Calculate ΔTm (use average of forward and reverse)
        if reference_seq:
            delta_tm_fwd = tm_target_fwd - tm_actual_fwd
            delta_tm_rev = tm_target_rev - tm_actual_rev
            delta_tm_avg = (delta_tm_fwd + delta_tm_rev) / 2.0
        else:
            delta_tm_fwd = delta_tm_rev = delta_tm_avg = 0.0

        # Analyze primer structure (GC, dimers, hairpins)
        structure_fwd = analyze_primer_structure(forward_seq)
        structure_rev = analyze_primer_structure(reverse_seq)

        # Analyze amplicon GC content for longer products (>1kb)
        amplicon_gc = None
        if amplicon_length > 1000 and self.genome_fasta and forward_pos is not None and reverse_pos is not None:
            try:
                record = next(SeqIO.parse(self.genome_fasta, 'fasta'))
                # Extract amplicon sequence between primers
                amplicon_seq = str(record.seq[forward_pos:reverse_pos])
                amplicon_gc = calculate_gc_content(amplicon_seq)
            except Exception as e:
                print(f"Warning: Could not calculate amplicon GC for {primer_name}: {e}", file=sys.stderr)

        # Perform precise alignment if genome FASTA available
        alignment_fwd = None
        alignment_rev = None

        if self.genome_fasta and forward_pos is not None and reverse_pos is not None:
            # Extract genomic sequences at binding sites
            # Use biological length (count [ABC] as 1 base)
            fwd_length = count_primer_length(forward_seq)
            rev_length = count_primer_length(reverse_seq)

            genomic_fwd = extract_genomic_sequence(
                self.genome_fasta, forward_pos, fwd_length, 'forward'
            )
            genomic_rev = extract_genomic_sequence(
                self.genome_fasta, reverse_pos, rev_length, 'reverse'
            )

            # Perform alignments if extraction succeeded
            if genomic_fwd:
                alignment_fwd = align_primer_to_genome(forward_seq, genomic_fwd)
            if genomic_rev:
                alignment_rev = align_primer_to_genome(reverse_seq, genomic_rev)

        # Analyze mismatch profiles
        # Use precise 3' mismatch counts from alignment if available
        # Note: Alignment counts may differ from primersearch because primersearch uses IUPAC ambiguity matching
        if alignment_fwd and 'mismatch_count' in alignment_fwd:
            actual_fwd_mm = alignment_fwd['mismatch_count']
            mm_profile_fwd = {
                'total_mismatches': actual_fwd_mm,
                'primer_length': count_primer_length(forward_seq),
                'estimated_3prime_mismatches': alignment_fwd['3prime_mismatch_count'],
                '3prime_region_size': 5,
                'mismatch_density': round(actual_fwd_mm / count_primer_length(forward_seq), 3) if count_primer_length(forward_seq) > 0 else 0,
                'precise': True
            }
        else:
            actual_fwd_mm = forward_mm
            mm_profile_fwd = self.analyze_mismatch_profile(forward_seq, forward_mm, 'forward')
            mm_profile_fwd['precise'] = False

        if alignment_rev and 'mismatch_count' in alignment_rev:
            actual_rev_mm = alignment_rev['mismatch_count']
            mm_profile_rev = {
                'total_mismatches': actual_rev_mm,
                'primer_length': count_primer_length(reverse_seq),
                'estimated_3prime_mismatches': alignment_rev['3prime_mismatch_count'],
                '3prime_region_size': 5,
                'mismatch_density': round(actual_rev_mm / count_primer_length(reverse_seq), 3) if count_primer_length(reverse_seq) > 0 else 0,
                'precise': True
            }
        else:
            actual_rev_mm = reverse_mm
            mm_profile_rev = self.analyze_mismatch_profile(reverse_seq, reverse_mm, 'reverse')
            mm_profile_rev['precise'] = False

        # Calculate amplification probability
        p_amp_fwd, class_fwd = self.calculate_amplification_probability(
            forward_seq, tm_actual_fwd, tm_target_fwd, delta_tm_fwd, mm_profile_fwd, 'forward', is_target)
        p_amp_rev, class_rev = self.calculate_amplification_probability(
            reverse_seq, tm_actual_rev, tm_target_rev, delta_tm_rev, mm_profile_rev, 'reverse', is_target)

        # Overall amplification probability (both primers must work)
        p_amp_overall = p_amp_fwd * p_amp_rev

        # Classification: intended targets take priority
        if class_fwd == 'target' and class_rev == 'target':
            classification = 'target'
        # For off-targets, classify based on OVERALL probability
        # This is more accurate than worst-case individual primer classification
        elif p_amp_overall >= 0.7:
            classification = 'high'
        elif p_amp_overall >= 0.3:
            classification = 'medium'
        else:
            classification = 'low'

        result = {
            'primer_name': primer_name,
            'genome': genome_name,
            'amplicon_length': amplicon_length,
            'is_intended_target': is_target,

            # Forward primer thermodynamics
            'forward_primer_seq': forward_seq,  # Original with IUPAC codes
            'forward_mismatches': actual_fwd_mm,  # Use alignment count if available
            'forward_mismatches_primersearch': forward_mm,  # Original primersearch count
            'forward_position': forward_pos,
            'tm_target_forward': tm_target_fwd,
            'tm_actual_forward': tm_actual_fwd,
            'delta_tm_forward': round(delta_tm_fwd, 2),
            'forward_3prime_mm_est': mm_profile_fwd['estimated_3prime_mismatches'],
            'forward_3prime_mm_precise': mm_profile_fwd.get('precise', False),

            # Reverse primer thermodynamics
            'reverse_primer_seq': reverse_seq,  # Original with IUPAC codes
            'reverse_mismatches': actual_rev_mm,  # Use alignment count if available
            'reverse_mismatches_primersearch': reverse_mm,  # Original primersearch count
            'reverse_position': reverse_pos,
            'tm_target_reverse': tm_target_rev,
            'tm_actual_reverse': tm_actual_rev,
            'delta_tm_reverse': round(delta_tm_rev, 2),
            'reverse_3prime_mm_est': mm_profile_rev['estimated_3prime_mismatches'],
            'reverse_3prime_mm_precise': mm_profile_rev.get('precise', False),

            # Overall metrics
            'delta_tm_average': round(delta_tm_avg, 2),
            'total_mismatches': actual_fwd_mm + actual_rev_mm,  # Use alignment counts
            'total_mismatches_primersearch': forward_mm + reverse_mm,  # Original primersearch count
            'p_amplification_forward': round(p_amp_fwd, 3),
            'p_amplification_reverse': round(p_amp_rev, 3),
            'p_amplification_overall': round(p_amp_overall, 3),
            'risk_classification': classification,

            # Primer structure analysis
            'forward_gc_content': structure_fwd['gc_content'],
            'forward_self_dimer_risk': structure_fwd['self_dimer_risk'],
            'forward_hairpin_risk': structure_fwd['hairpin_risk'],
            'reverse_gc_content': structure_rev['gc_content'],
            'reverse_self_dimer_risk': structure_rev['self_dimer_risk'],
            'reverse_hairpin_risk': structure_rev['hairpin_risk'],

            # Amplicon metrics
            'amplicon_gc_content': amplicon_gc if amplicon_gc is not None else 'N/A',

            # PCR conditions used
            'annealing_temp': self.annealing_temp,
            'buffer_Na_mM': self.pcr_conditions.get('Na', 50.0),
            'buffer_Mg_mM': self.pcr_conditions.get('Mg', 2.0)
        }

        # Add alignment data if available
        if alignment_fwd:
            result['forward_alignment'] = {
                'alignment_text': alignment_fwd.get('alignment_text', ''),
                'mismatch_positions': alignment_fwd.get('mismatch_positions', []),
                'mismatch_details': alignment_fwd.get('mismatch_details', []),
                'mismatch_count': alignment_fwd.get('mismatch_count', 0),
                '3prime_mismatch_count': alignment_fwd.get('3prime_mismatch_count', 0),
                'expanded_primer': alignment_fwd.get('expanded_primer', '')
            }

        if alignment_rev:
            result['reverse_alignment'] = {
                'alignment_text': alignment_rev.get('alignment_text', ''),
                'mismatch_positions': alignment_rev.get('mismatch_positions', []),
                'mismatch_details': alignment_rev.get('mismatch_details', []),
                'mismatch_count': alignment_rev.get('mismatch_count', 0),
                '3prime_mismatch_count': alignment_rev.get('3prime_mismatch_count', 0),
                'expanded_primer': alignment_rev.get('expanded_primer', '')
            }

        return result


def parse_primersearch_output(primersearch_file: Path) -> List[Dict]:
    """
    Parse EMBOSS primersearch output to extract amplicon data

    Returns:
        List of dicts with amplicon information
    """
    amplicons = []
    current_primer = None
    current_amplicon = {}

    with open(primersearch_file) as f:
        for line in f:
            line = line.rstrip()

            # Primer name line
            if line.startswith('Primer name'):
                match = re.search(r'Primer name (.+)', line)
                if match:
                    current_primer = match.group(1)

            # Amplimer header
            elif line.startswith('Amplimer'):
                if current_amplicon:
                    amplicons.append(current_amplicon)
                current_amplicon = {'primer_name': current_primer}

            # Sequence line
            elif line.strip().startswith('Sequence:'):
                match = re.search(r'Sequence:\s+(\S+)', line)
                if match:
                    current_amplicon['genome'] = match.group(1)

            # Primer hit lines - extract sequence and mismatches
            elif 'hits' in line and 'at' in line:
                # Example: AGAGTTTGATCMTGGCTCAG hits forward strand at 9292 with 0 mismatches
                # Position can be plain (9292) or in brackets ([432702]) for reverse strand
                match = re.search(r'(\S+)\s+hits\s+(\w+)\s+strand\s+at\s+\[?(\d+)\]?\s+with\s+(\d+)\s+mismatch', line)
                if match:
                    primer_seq = match.group(1)
                    strand = match.group(2)
                    position = int(match.group(3))
                    mismatches = int(match.group(4))

                    if strand == 'forward':
                        current_amplicon['forward_seq'] = primer_seq
                        current_amplicon['forward_pos'] = position
                        current_amplicon['forward_mm'] = mismatches
                    else:
                        current_amplicon['reverse_seq'] = primer_seq
                        current_amplicon['reverse_pos'] = position
                        current_amplicon['reverse_mm'] = mismatches

            # Amplicon length
            elif line.strip().startswith('Amplimer length:'):
                match = re.search(r'Amplimer length:\s+(\d+)', line)
                if match:
                    current_amplicon['amplicon_length'] = int(match.group(1))

        # Add last amplicon
        if current_amplicon:
            amplicons.append(current_amplicon)

    return amplicons


def extract_genomic_sequence(genome_fasta: Path, position: int, length: int, strand: str) -> Optional[str]:
    """
    Extract genomic sequence at primer binding site

    Args:
        genome_fasta: Path to genome FASTA file
        position: Binding position (1-indexed, as reported by primersearch)
        length: Length to extract (primer length after cleaning IUPAC codes)
        strand: 'forward' or 'reverse'

    Returns:
        Genomic sequence at binding site (in primer orientation), or None if failed
    """
    try:
        # Load genome sequence (assume single-record FASTA for now)
        record = next(SeqIO.parse(genome_fasta, 'fasta'))

        # Convert to 0-indexed
        pos_0idx = position - 1

        if strand == 'forward':
            # Extract forward sequence
            binding_seq = str(record.seq[pos_0idx:pos_0idx + length])
        else:
            # For reverse strand, primersearch reports the end position
            # We need to extract the sequence and reverse complement it
            binding_seq = str(record.seq[pos_0idx - length + 1:pos_0idx + 1].reverse_complement())

        return binding_seq

    except Exception as e:
        print(f"Warning: Failed to extract genomic sequence at position {position}: {e}", file=sys.stderr)
        return None


def count_primer_length(primer_seq: str) -> int:
    """
    Count the biological length of a primer sequence.
    IUPAC codes like [ABC] count as 1 base, ? counts as 1 base.

    Args:
        primer_seq: Primer sequence with potential IUPAC codes

    Returns:
        Biological length in bases
    """
    # Count [ABC] as 1, ? as 1
    length = 0
    i = 0
    while i < len(primer_seq):
        if primer_seq[i] == '[':
            # Skip until ]
            while i < len(primer_seq) and primer_seq[i] != ']':
                i += 1
            length += 1
            i += 1
        else:
            length += 1
            i += 1
    return length


def expand_primer_to_match_genome(primer_seq: str, genomic_seq: str) -> str:
    """
    Expand degenerate primer by selecting variants that match the genomic sequence.

    For each IUPAC ambiguity code in the primer, choose the variant that matches
    the corresponding position in the genomic sequence. This ensures that when
    primersearch reports 0 mismatches, the alignment will also show 0 mismatches.

    Args:
        primer_seq: Primer sequence with IUPAC codes (e.g., "AG[CAM]GTT[TCY]")
        genomic_seq: Genomic sequence at binding site

    Returns:
        Expanded primer sequence using matching variants (e.g., "AGCGTTT")
    """
    expanded = []
    primer_bio_pos = 0  # Position in biological primer (counting [ABC] as 1)
    i = 0  # Position in primer string

    while i < len(primer_seq):
        if primer_seq[i] == '[':
            # Find the closing bracket
            close_idx = primer_seq.index(']', i)
            variants = primer_seq[i+1:close_idx]

            # Check which variant matches genome at this position
            if primer_bio_pos < len(genomic_seq):
                genome_base = genomic_seq[primer_bio_pos].upper()
                if genome_base in variants:
                    # Use the matching variant
                    expanded.append(genome_base)
                else:
                    # No match - take first variant (shouldn't happen for intended targets)
                    expanded.append(variants[0])
            else:
                # Genomic sequence too short - take first variant
                expanded.append(variants[0])

            primer_bio_pos += 1
            i = close_idx + 1

        elif primer_seq[i] == '?':
            # Wildcard - use genome base if available
            if primer_bio_pos < len(genomic_seq):
                expanded.append(genomic_seq[primer_bio_pos])
            else:
                expanded.append('N')
            primer_bio_pos += 1
            i += 1

        else:
            # Regular base
            expanded.append(primer_seq[i])
            primer_bio_pos += 1
            i += 1

    return ''.join(expanded)


def align_primer_to_genome(primer_seq: str, genomic_seq: str) -> Dict:
    """
    Perform pairwise alignment between primer and genomic sequence

    Uses BioPython's PairwiseAligner to create detailed alignment showing
    sequence identity (not complementary base pairing).

    Args:
        primer_seq: Primer sequence (may contain IUPAC ambiguity codes)
        genomic_seq: Actual genomic sequence at binding site

    Returns:
        Dict with alignment string, mismatch positions, and mismatch counts
    """
    # Expand IUPAC ambiguity codes to match the genomic sequence
    # This ensures intended targets show 0 mismatches when primersearch reports 0
    expanded_primer = expand_primer_to_match_genome(primer_seq, genomic_seq)

    # Also handle wildcards
    clean_primer = re.sub(r'\?', 'N', expanded_primer)

    # Handle cases where sequences don't match in length
    if len(genomic_seq) == 0:
        return {
            'alignment_text': '',
            'mismatch_positions': [],
            'mismatch_count': 0,
            '3prime_mismatch_count': 0,
            'alignment_score': 0
        }

    # Create aligner - use global mode to align entire sequences
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    try:
        # Perform alignment
        alignments = aligner.align(clean_primer, genomic_seq)
        if len(alignments) == 0:
            raise ValueError("No alignment found")

        best_alignment = alignments[0]

        # Extract alignment strings
        primer_aligned = str(best_alignment[0])
        genome_aligned = str(best_alignment[1])

        # Analyze mismatches
        mismatches = []
        mismatch_positions = []

        # Track position in ungapped primer for accurate 3' distance calculation
        ungapped_primer_pos = 0

        for i, (p, g) in enumerate(zip(primer_aligned, genome_aligned)):
            # Count ungapped positions in primer
            if p != '-':
                current_primer_pos = ungapped_primer_pos
                ungapped_primer_pos += 1

            if p != g and p != '-' and g != '-':
                # Distance from 3' end in the ungapped primer
                distance_from_3prime = len(clean_primer) - current_primer_pos - 1

                mismatches.append({
                    'position': i,  # Position in aligned sequence
                    'primer_base': p,
                    'genome_base': g,
                    'distance_from_3prime': distance_from_3prime
                })
                mismatch_positions.append(i)

        # Count 3' mismatches (last 5 bases, excluding gaps)
        three_prime_threshold = 5
        # Count mismatches in last 5 non-gap positions from 3' end
        non_gap_positions = [m for m in mismatches if m['distance_from_3prime'] >= 0]
        three_prime_mm = sum(1 for m in non_gap_positions if m['distance_from_3prime'] < three_prime_threshold)

        # Count gaps
        gap_count = sum(1 for p, g in zip(primer_aligned, genome_aligned) if p == '-' or g == '-')

        # Generate alignment visualization text
        match_string = ''
        for p, g in zip(primer_aligned, genome_aligned):
            if p == g:
                match_string += '|'
            elif p == '-' or g == '-':
                match_string += ' '
            else:
                match_string += '.'

        alignment_text = f"{primer_aligned}\n{match_string}\n{genome_aligned}"

        # Add warning if there are many gaps (suggests extraction issue)
        if gap_count > 3:
            alignment_text += f"\nWarning: {gap_count} gaps in alignment - position extraction may be imprecise"

        return {
            'alignment_text': alignment_text,
            'primer_aligned': primer_aligned,
            'genome_aligned': genome_aligned,
            'match_string': match_string,
            'mismatch_positions': mismatch_positions,
            'mismatch_details': mismatches,
            'mismatch_count': len(mismatches),
            '3prime_mismatch_count': three_prime_mm,
            'alignment_score': float(best_alignment.score),
            'expanded_primer': expanded_primer  # Show which variant matched
        }

    except Exception as e:
        print(f"Warning: Alignment failed for {primer_seq[:20]}...: {e}", file=sys.stderr)
        return {
            'alignment_text': '',
            'mismatch_positions': [],
            'mismatch_count': 0,
            '3prime_mismatch_count': 0,
            'alignment_score': 0
        }


def calculate_gc_content(seq: str) -> float:
    """
    Calculate GC content as percentage

    Args:
        seq: DNA sequence

    Returns:
        GC content as percentage (0-100)
    """
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    total = len(seq_upper)
    return round((gc_count / total * 100), 2) if total > 0 else 0.0


def check_self_dimer(seq: str) -> Dict:
    """
    Check for potential self-dimer formation using simple complementarity

    This is a simplified check - looks for reverse-complementary regions
    that could form stable dimers.

    Args:
        seq: Primer sequence

    Returns:
        Dict with dimer score and details
    """
    # Create reverse complement
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    rev_comp = ''.join(complement_map.get(b, 'N') for b in reversed(seq.upper()))

    # Look for complementary regions at 3' end (most critical for extension)
    max_match = 0
    for i in range(len(seq) - 2):  # Check last 3+ bases
        matches = 0
        for j in range(min(len(seq) - i, 6)):  # Check up to 6 bases
            if i + j < len(seq) and j < len(rev_comp):
                if seq[-(i+j+1)] == rev_comp[j]:
                    matches += 1
                else:
                    break
        max_match = max(max_match, matches)

    # Score based on consecutive matches at 3' end
    # >4 consecutive matches at 3' end = high dimer risk
    dimer_risk = 'high' if max_match >= 4 else ('medium' if max_match == 3 else 'low')

    return {
        'max_3prime_complementarity': max_match,
        'dimer_risk': dimer_risk
    }


def check_hairpin(seq: str) -> Dict:
    """
    Check for potential hairpin formation using palindrome detection

    Simplified hairpin check - looks for inverted repeats that could
    form stable secondary structures.

    Args:
        seq: Primer sequence

    Returns:
        Dict with hairpin score and details
    """
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    max_stem_length = 0
    hairpin_positions = []

    # Check for palindromic regions (potential hairpin stems)
    for i in range(len(seq) - 3):  # Minimum stem = 3bp
        for j in range(i + 4, min(i + 12, len(seq))):  # Loop 4-12 bases
            stem_len = 0
            # Check for reverse complement after the loop
            for k in range(min(len(seq) - j, j - i - 3)):
                if seq[i + k].upper() in complement_map:
                    comp_base = complement_map[seq[i + k].upper()]
                    if j + k < len(seq) and seq[j + k].upper() == comp_base:
                        stem_len += 1
                    else:
                        break
            if stem_len >= 3:
                max_stem_length = max(max_stem_length, stem_len)
                hairpin_positions.append((i, j, stem_len))

    # Score based on stem length
    # >5bp stem = high hairpin risk (very stable)
    # 4-5bp stem = medium risk
    # <4bp stem = low risk
    hairpin_risk = 'high' if max_stem_length >= 5 else ('medium' if max_stem_length >= 4 else 'low')

    return {
        'max_stem_length': max_stem_length,
        'hairpin_risk': hairpin_risk,
        'hairpin_count': len(hairpin_positions)
    }


def analyze_primer_structure(seq: str) -> Dict:
    """
    Comprehensive primer structure analysis

    Args:
        seq: Primer sequence

    Returns:
        Dict with GC content and secondary structure predictions
    """
    gc_content = calculate_gc_content(seq)
    self_dimer = check_self_dimer(seq)
    hairpin = check_hairpin(seq)

    return {
        'gc_content': gc_content,
        'self_dimer_score': self_dimer['max_3prime_complementarity'],
        'self_dimer_risk': self_dimer['dimer_risk'],
        'hairpin_stem_length': hairpin['max_stem_length'],
        'hairpin_risk': hairpin['hairpin_risk']
    }


def load_pcr_conditions(master_mix_file: Path, master_mix_name: str) -> Dict:
    """Load PCR conditions from master mix JSON file"""
    with open(master_mix_file) as f:
        master_mixes = json.load(f)

    if master_mix_name not in master_mixes:
        print(f"Error: Master mix '{master_mix_name}' not found in {master_mix_file}", file=sys.stderr)
        print(f"Available mixes: {', '.join([k for k in master_mixes.keys() if not k.startswith('_')])}",
              file=sys.stderr)
        sys.exit(1)

    mix = master_mixes[master_mix_name]
    return mix['buffer_composition']


def load_primers(primer_file: Path) -> Dict[str, Tuple[str, str]]:
    """
    Load primers from EMBOSS primersearch format

    Returns:
        Dict mapping primer_name -> (forward_seq, reverse_seq)
    """
    primers = {}
    with open(primer_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 3:
                primer_name = parts[0]
                forward_seq = parts[1]
                reverse_seq = parts[2]
                primers[primer_name] = (forward_seq, reverse_seq)

    return primers


def load_references(reference_file: Path) -> Dict[str, str]:
    """
    Load reference target sequences from FASTA

    Returns:
        Dict mapping primer_name (from FASTA header) -> sequence
    """
    references = {}

    if not reference_file.exists():
        print(f"Warning: Reference file {reference_file} not found", file=sys.stderr)
        return references

    for record in SeqIO.parse(reference_file, 'fasta'):
        # Extract primer name from header (first word after >)
        primer_name = record.id.split('_reference')[0]
        references[primer_name] = str(record.seq)

    return references


def main():
    parser = argparse.ArgumentParser(
        description='Thermodynamic analysis of primer binding sites',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('primersearch_file', type=Path,
                       help='Input primersearch output file')
    parser.add_argument('primers_file', type=Path,
                       help='Primers file (EMBOSS format)')
    parser.add_argument('--genome-fasta', type=Path,
                       help='Genome FASTA file for extracting binding sequences and alignments')
    parser.add_argument('--reference', type=Path,
                       help='Reference target sequences (FASTA)')
    parser.add_argument('--genome-name', type=str,
                       help='Original genome query name (overrides name from primersearch output)')
    parser.add_argument('--master-mix', default='DreamTaq',
                       help='Master mix name from master_mixes.json (default: DreamTaq)')
    parser.add_argument('--master-mix-file', type=Path,
                       default=Path(__file__).parent.parent / 'pcr_conditions/master_mixes.json',
                       help='Path to master_mixes.json')
    parser.add_argument('--annealing-temp', type=float, default=60.0,
                       help='Annealing temperature in °C (default: 60)')
    parser.add_argument('--delta-tm-high', type=float, default=5.0,
                       help='ΔTm threshold for high amplification risk (default: 5°C)')
    parser.add_argument('--delta-tm-medium', type=float, default=10.0,
                       help='ΔTm threshold for medium amplification risk (default: 10°C)')
    parser.add_argument('-o', '--output', type=Path,
                       help='Output JSON file (default: stdout)')

    args = parser.parse_args()

    # Load inputs
    print(f"Loading PCR conditions: {args.master_mix}", file=sys.stderr)
    pcr_conditions = load_pcr_conditions(args.master_mix_file, args.master_mix)

    print(f"Loading primers from {args.primers_file}", file=sys.stderr)
    primers = load_primers(args.primers_file)

    references = {}
    if args.reference:
        print(f"Loading reference sequences from {args.reference}", file=sys.stderr)
        references = load_references(args.reference)

    print(f"Parsing primersearch results from {args.primersearch_file}", file=sys.stderr)
    amplicons = parse_primersearch_output(args.primersearch_file)

    # Initialize analyzer
    analyzer = ThermoAnalyzer(
        pcr_conditions=pcr_conditions,
        primers=primers,
        references=references,
        genome_fasta=args.genome_fasta,
        annealing_temp=args.annealing_temp,
        delta_tm_thresholds=(args.delta_tm_high, args.delta_tm_medium)
    )

    # Analyze each amplicon
    print(f"Analyzing {len(amplicons)} amplicons...", file=sys.stderr)
    if args.genome_fasta:
        print(f"Using genome FASTA for precise alignment analysis: {args.genome_fasta}", file=sys.stderr)
    results = []

    # Use provided genome name or extract from primersearch output
    genome_display_name = args.genome_name if args.genome_name else None

    for amp in amplicons:
        if not all(k in amp for k in ['forward_seq', 'reverse_seq', 'forward_mm', 'reverse_mm']):
            continue  # Skip incomplete amplicons

        # Use provided genome name or fall back to name from primersearch output
        genome_name = genome_display_name if genome_display_name else amp.get('genome', 'Unknown')

        analysis = analyzer.analyze_amplicon(
            primer_name=amp['primer_name'],
            genome_name=genome_name,
            forward_seq=amp['forward_seq'],
            reverse_seq=amp['reverse_seq'],
            forward_mm=amp['forward_mm'],
            reverse_mm=amp['reverse_mm'],
            amplicon_length=amp.get('amplicon_length', 0),
            forward_pos=amp.get('forward_pos'),
            reverse_pos=amp.get('reverse_pos')
        )
        results.append(analysis)

    # Output results
    output_data = {
        'metadata': {
            'primersearch_file': str(args.primersearch_file),
            'master_mix': args.master_mix,
            'annealing_temp': args.annealing_temp,
            'pcr_conditions': pcr_conditions,
            'delta_tm_thresholds': {
                'high_risk': args.delta_tm_high,
                'medium_risk': args.delta_tm_medium
            },
            'total_amplicons': len(results)
        },
        'results': results
    }

    if args.output:
        with open(args.output, 'w') as f:
            json.dump(output_data, f, indent=2)
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(json.dumps(output_data, indent=2))

    # Summary statistics
    high_risk = sum(1 for r in results if r['risk_classification'] == 'high')
    medium_risk = sum(1 for r in results if r['risk_classification'] == 'medium')
    low_risk = sum(1 for r in results if r['risk_classification'] == 'low')

    print(f"\nAnalysis complete:", file=sys.stderr)
    print(f"  High risk amplifications: {high_risk}", file=sys.stderr)
    print(f"  Medium risk amplifications: {medium_risk}", file=sys.stderr)
    print(f"  Low risk amplifications: {low_risk}", file=sys.stderr)


if __name__ == '__main__':
    main()
