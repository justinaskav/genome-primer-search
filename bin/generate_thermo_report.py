#!/usr/bin/env python3
"""
Generate Thermodynamic Analysis Reports

Creates comprehensive reports from thermodynamic analysis JSON:
- Detailed per-amplicon TSV
- Off-target summary
- Primer specificity metrics
- Interactive HTML report
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import List, Dict


def is_target(r: Dict) -> bool:
    """Hit was matched against a declared reference and accepted as intended."""
    return r.get('is_intended_target') is True


def is_off_target(r: Dict) -> bool:
    """Hit was matched against a declared reference and rejected."""
    return r.get('is_intended_target') is False


def is_unclassified(r: Dict) -> bool:
    """No reference was declared for this primer — we make no claim."""
    return r.get('is_intended_target') is None


def target_label(r: Dict) -> str:
    if is_target(r):
        return 'Yes'
    if is_off_target(r):
        return 'No'
    return 'Undeclared'


def _intent_banner_html(references_provided: bool,
                        references_per_primer: Dict[str, bool]) -> str:
    """
    Banner at the top of the HTML report stating whether the user declared
    intent (reference FASTA) and which primers it covers. Without this,
    readers can't tell whether the absence of an "Off-Targets" section means
    "no off-targets found" or "intent was never declared."
    """
    if not references_provided:
        return (
            '<div style="background:#fff8e1;border-left:4px solid #f9a825;'
            'padding:15px;margin:20px 0;">'
            '<strong>Intent declaration: none.</strong> '
            'No reference FASTA was provided, so this report cannot classify '
            'hits as intended targets vs. off-targets. Tm, ΔTm, P(amplification) '
            'and mismatch profiles below are computed from primer + binding-site '
            'sequence alone and are still valid — interpret them comparatively. '
            'To enable target/off-target classification, pass '
            '<code>--thermo_references &lt;file.fasta&gt;</code> with one entry '
            'per primer (header <code>&gt;primer_name_reference</code>).'
            '</div>'
        )
    covered = sum(1 for v in references_per_primer.values() if v)
    total = len(references_per_primer)
    missing = [p for p, v in references_per_primer.items() if not v]
    missing_note = ''
    if missing:
        missing_note = (
            ' Hits from primers without a reference will appear as '
            '<em>Undeclared</em> and are excluded from off-target counts.'
        )
        if len(missing) <= 6:
            missing_note += f' Missing: {", ".join(sorted(missing))}.'
    return (
        '<div style="background:#e8f5e9;border-left:4px solid #43a047;'
        'padding:15px;margin:20px 0;">'
        f'<strong>Intent declaration:</strong> reference provided for '
        f'{covered} of {total} primers. Intended targets are hits whose '
        f'binding region matches the declared reference at &ge; 95 % identity.'
        f'{missing_note}'
        '</div>'
    )


def generate_detailed_tsv(results: List[Dict], output_file: Path):
    """Generate detailed per-amplicon TSV report"""

    headers = [
        'Primer', 'Genome', 'Forward_Position', 'Reverse_Position', 'Amplicon_Length',
        'Is_Target', 'Risk_Classification', 'P_Amplification',
        'Forward_Tm_Target', 'Forward_Tm_Actual', 'Forward_Delta_Tm', 'Forward_Mismatches',
        'Reverse_Tm_Target', 'Reverse_Tm_Actual', 'Reverse_Delta_Tm', 'Reverse_Mismatches',
        'Delta_Tm_Average', 'Total_Mismatches',
        'Annealing_Temp', 'Buffer_Na_mM', 'Buffer_Mg_mM'
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

        for r in results:
            row = [
                r['primer_name'],
                r['genome'],
                r.get('forward_position', 'N/A'),
                r.get('reverse_position', 'N/A'),
                r['amplicon_length'],
                target_label(r),
                r['risk_classification'],
                r['p_amplification_overall'],
                r['tm_target_forward'],
                r['tm_actual_forward'],
                r['delta_tm_forward'],
                r['forward_mismatches'],
                r['tm_target_reverse'],
                r['tm_actual_reverse'],
                r['delta_tm_reverse'],
                r['reverse_mismatches'],
                r['delta_tm_average'],
                r['total_mismatches'],
                r['annealing_temp'],
                r['buffer_Na_mM'],
                r['buffer_Mg_mM']
            ]
            f.write('\t'.join(map(str, row)) + '\n')


def generate_offtarget_summary(results: List[Dict], output_file: Path):
    """Generate summary of high-risk off-target amplifications"""

    headers = [
        'Primer', 'Genome', 'Forward_Position', 'Reverse_Position',
        'Amplicon_Length', 'Risk_Classification', 'P_Amplification',
        'Delta_Tm_Average', 'Total_Mismatches',
        'Forward_3prime_MM', 'Reverse_3prime_MM',
        'Recommendation'
    ]

    # Off-target list only contains hits that were *declared* off-target
    # (reference exists for the primer and the binding region did not match).
    # Undeclared hits (no reference for this primer) are excluded — calling
    # them "off-target" with no intent declared would be a lie.
    offtargets = [r for r in results
                  if is_off_target(r)
                  and r['risk_classification'] in ['high', 'medium']]

    # Sort by risk (high first) then by P_amplification (descending)
    offtargets.sort(key=lambda x: (
        0 if x['risk_classification'] == 'high' else 1,
        -x['p_amplification_overall']
    ))

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

        for r in offtargets:
            # Generate recommendation
            if r['risk_classification'] == 'high':
                recommendation = 'REVIEW REQUIRED - High amplification probability'
            elif r['delta_tm_average'] < 5:
                recommendation = 'CAUTION - Low ΔTm, potential off-target'
            elif r['total_mismatches'] == 0:
                recommendation = 'WARNING - Perfect match to off-target'
            else:
                recommendation = 'Monitor - Medium risk off-target'

            row = [
                r['primer_name'],
                r['genome'],
                r.get('forward_position', 'N/A'),
                r.get('reverse_position', 'N/A'),
                r['amplicon_length'],
                r['risk_classification'].upper(),
                r['p_amplification_overall'],
                r['delta_tm_average'],
                r['total_mismatches'],
                r['forward_3prime_mm_est'],
                r['reverse_3prime_mm_est'],
                recommendation
            ]
            f.write('\t'.join(map(str, row)) + '\n')


def generate_primer_specificity(results: List[Dict], output_file: Path):
    """
    Generate per-primer specificity metrics.

    Specificity is only meaningful when the user declared intent (reference
    FASTA). Without that, the score collapses to nonsense — 0 targets and 0
    off-targets yields a vacuous "100% Excellent" rating. We detect the
    no-intent case by checking whether any result has a non-null
    is_intended_target and, if not, write a header-only TSV so downstream
    consumers see the same intent state the HTML conveys.
    """
    intent_declared = any(r.get('is_intended_target') is not None for r in results)
    if not intent_declared:
        with open(output_file, 'w') as f:
            f.write(
                'Primer\tTotal_Hits\tTarget_Hits\tOffTarget_Hits\tHigh_Risk_OffTargets\t'
                'Medium_Risk_OffTargets\tLow_Risk_OffTargets\tGenomes_Hit\t'
                'Mean_OffTarget_Delta_Tm\tMean_OffTarget_P_Amp\tSpecificity_Score\tSpecificity_Rating\n'
            )
        return

    # Aggregate by primer
    primer_stats = defaultdict(lambda: {
        'total_hits': 0,
        'target_hits': 0,
        'offtarget_hits': 0,
        'high_risk_offtargets': 0,
        'medium_risk_offtargets': 0,
        'low_risk_offtargets': 0,
        'genomes_hit': set(),
        'mean_delta_tm': [],
        'mean_p_amp_offtarget': []
    })

    for r in results:
        primer = r['primer_name']
        stats = primer_stats[primer]

        stats['total_hits'] += 1
        stats['genomes_hit'].add(r['genome'])

        if is_target(r):
            stats['target_hits'] += 1
        elif is_off_target(r):
            stats['offtarget_hits'] += 1
            stats['mean_delta_tm'].append(r['delta_tm_average'])
            stats['mean_p_amp_offtarget'].append(r['p_amplification_overall'])

            if r['risk_classification'] == 'high':
                stats['high_risk_offtargets'] += 1
            elif r['risk_classification'] == 'medium':
                stats['medium_risk_offtargets'] += 1
            else:
                stats['low_risk_offtargets'] += 1
        # Undeclared hits (no per-primer reference) are counted in total_hits
        # only — they don't pollute the target / off-target tallies.

    headers = [
        'Primer', 'Total_Hits', 'Target_Hits', 'OffTarget_Hits',
        'High_Risk_OffTargets', 'Medium_Risk_OffTargets', 'Low_Risk_OffTargets',
        'Genomes_Hit', 'Mean_OffTarget_Delta_Tm', 'Mean_OffTarget_P_Amp',
        'Specificity_Score', 'Specificity_Rating'
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

        for primer, stats in sorted(primer_stats.items()):
            # Calculate mean ΔTm for off-targets
            if stats['mean_delta_tm']:
                mean_delta_tm = round(sum(stats['mean_delta_tm']) / len(stats['mean_delta_tm']), 2)
            else:
                mean_delta_tm = 'N/A'

            # Calculate mean P_amp for off-targets
            if stats['mean_p_amp_offtarget']:
                mean_p_amp = round(sum(stats['mean_p_amp_offtarget']) / len(stats['mean_p_amp_offtarget']), 3)
            else:
                mean_p_amp = 'N/A'

            # Calculate specificity score (0-100)
            # Higher is better: Weighted by off-target risk levels
            # High risk = 1.0 weight, Medium = 0.5 weight, Low = 0.1 weight
            if stats['total_hits'] > 0:
                weighted_offtargets = (
                    stats['high_risk_offtargets'] * 1.0 +
                    stats['medium_risk_offtargets'] * 0.5 +
                    stats['low_risk_offtargets'] * 0.1
                )
                specificity_score = round(
                    100 * (1 - (weighted_offtargets / stats['total_hits'])),
                    1
                )
            else:
                specificity_score = 0

            # Rating
            if specificity_score >= 90:
                rating = 'Excellent'
            elif specificity_score >= 75:
                rating = 'Good'
            elif specificity_score >= 50:
                rating = 'Fair'
            else:
                rating = 'Poor'

            row = [
                primer,
                stats['total_hits'],
                stats['target_hits'],
                stats['offtarget_hits'],
                stats['high_risk_offtargets'],
                stats['medium_risk_offtargets'],
                stats['low_risk_offtargets'],
                len(stats['genomes_hit']),
                mean_delta_tm,
                mean_p_amp,
                specificity_score,
                rating
            ]
            f.write('\t'.join(map(str, row)) + '\n')


def format_position_range(fwd_pos, rev_pos, amp_len=None):
    """
    Format primer position range with calculated amplicon coordinates.

    Args:
        fwd_pos: Forward primer position (start)
        rev_pos: Reverse primer position
        amp_len: Amplicon length (optional, for calculating actual range)

    Returns:
        Formatted position string showing amplicon start..end
    """
    if not fwd_pos:
        return "N/A"

    if amp_len:
        # Calculate actual amplicon end from forward position + length
        amp_start = fwd_pos
        amp_end = fwd_pos + amp_len - 1

        if fwd_pos > rev_pos:
            # Circular genome: amplicon wraps around origin
            return f"{amp_start:,}..{amp_end:,} (⟳)"
        else:
            # Linear region
            return f"{amp_start:,}..{amp_end:,}"
    else:
        # Fallback: use primer positions
        if not rev_pos:
            return f"{fwd_pos:,}.."
        elif fwd_pos > rev_pos:
            return f"{fwd_pos:,}→{rev_pos:,} (⟳)"
        else:
            return f"{fwd_pos:,}..{rev_pos:,}"


def format_alignment_html(result: Dict, row_id: str) -> str:
    """
    Format alignment data as HTML for expandable row

    Args:
        result: Amplicon result dict with alignment data
        row_id: Unique row identifier

    Returns:
        HTML string with alignment visualization
    """
    if 'forward_alignment' not in result and 'reverse_alignment' not in result:
        return ""

    html = f'<div class="alignment-container">'

    # Show genomic coordinates
    fwd_pos = result.get('forward_position')
    rev_pos = result.get('reverse_position')
    genome = result.get('genome', 'Unknown')
    amp_len = result.get('amplicon_length', 0)

    if fwd_pos and rev_pos and amp_len:
        html += f'<div style="margin-bottom: 10px; font-size: 12px; color: #1976d2; background: #e3f2fd; padding: 5px; border-left: 3px solid #2196F3;">'

        # Calculate amplicon range from forward position + length
        amp_start = fwd_pos
        amp_end = fwd_pos + amp_len - 1

        # Detect circular genome (forward > reverse)
        if fwd_pos > rev_pos:
            # Circular: amplicon wraps around origin
            html += f'<strong>Amplicon Region:</strong> {genome}:{amp_start:,}..{amp_end:,} '
            html += f'<span style="color: #d32f2f; font-weight: bold;">(circular, {amp_len}bp)</span><br>'
            html += f'<span style="font-size: 11px;">Wraps around genome origin | Fwd primer@{fwd_pos:,}, Rev primer@{rev_pos:,}</span>'
        else:
            # Linear: normal amplicon
            html += f'<strong>Amplicon Region:</strong> {genome}:{amp_start:,}..{amp_end:,} ({amp_len}bp)<br>'
            html += f'<span style="font-size: 11px;">Fwd primer@{fwd_pos:,}, Rev primer@{rev_pos:,}</span>'

        html += '</div>'

    # Show original primers with IUPAC codes
    fwd_primer_orig = result.get('forward_primer_seq', '')
    rev_primer_orig = result.get('reverse_primer_seq', '')

    html += f'<div style="margin-bottom: 10px; font-size: 12px; color: #666;">'
    html += f'<strong>Original Primers (degenerate):</strong><br>'
    html += f'Forward: <code>{fwd_primer_orig}</code><br>'
    html += f'Reverse: <code>{rev_primer_orig}</code>'
    html += '</div>'

    # Show matched variants (what was actually aligned)
    fwd_expanded = ''
    rev_expanded = ''
    if 'forward_alignment' in result:
        fwd_expanded = result['forward_alignment'].get('expanded_primer', '')
    if 'reverse_alignment' in result:
        rev_expanded = result['reverse_alignment'].get('expanded_primer', '')

    if fwd_expanded or rev_expanded:
        html += f'<div style="margin-bottom: 10px; font-size: 12px; color: #2c5f2d; background: #e8f5e9; padding: 5px; border-left: 3px solid #4caf50;">'
        html += f'<strong>Matched Variant (used for alignment):</strong><br>'
        if fwd_expanded:
            html += f'Forward: <code>{fwd_expanded}</code><br>'
        if rev_expanded:
            html += f'Reverse: <code>{rev_expanded}</code>'
        html += '</div>'

    # Forward primer alignment
    if 'forward_alignment' in result and result['forward_alignment'].get('alignment_text'):
        fwd_align = result['forward_alignment']
        # Use the precise mismatch counts from alignment
        fwd_mm = fwd_align.get('mismatch_count', 0)
        fwd_3prime_mm = fwd_align.get('3prime_mismatch_count', 0)
        fwd_mm_ps = result.get('forward_mismatches_primersearch', fwd_mm)

        # Color code based on mismatches
        if fwd_mm == 0:
            mm_class = 'match-good'
        elif fwd_mm <= 2:
            mm_class = 'match-warn'
        else:
            mm_class = 'match-bad'

        # Always show both counts for transparency and consistency
        mm_display = f'{fwd_mm} variant mismatch{"es" if fwd_mm != 1 else ""} | {fwd_mm_ps} degenerate mismatch{"es" if fwd_mm_ps != 1 else ""}'

        html += f'''
        <div class="alignment-header {mm_class}">
            Forward: {mm_display} ({fwd_3prime_mm} at 3' end)
        </div>
        <div class="alignment-text">{fwd_align['alignment_text']}</div>
        '''

    # Reverse primer alignment
    if 'reverse_alignment' in result and result['reverse_alignment'].get('alignment_text'):
        rev_align = result['reverse_alignment']
        # Use the precise mismatch counts from alignment
        rev_mm = rev_align.get('mismatch_count', 0)
        rev_3prime_mm = rev_align.get('3prime_mismatch_count', 0)
        rev_mm_ps = result.get('reverse_mismatches_primersearch', rev_mm)

        # Color code based on mismatches
        if rev_mm == 0:
            mm_class = 'match-good'
        elif rev_mm <= 2:
            mm_class = 'match-warn'
        else:
            mm_class = 'match-bad'

        # Always show both counts for transparency and consistency
        mm_display = f'{rev_mm} variant mismatch{"es" if rev_mm != 1 else ""} | {rev_mm_ps} degenerate mismatch{"es" if rev_mm_ps != 1 else ""}'

        html += f'''
        <div class="alignment-header {mm_class}" style="margin-top: 15px;">
            Reverse: {mm_display} ({rev_3prime_mm} at 3' end)
        </div>
        <div class="alignment-text">{rev_align['alignment_text']}</div>
        '''

    html += '</div>'
    return html


def generate_html_report(results: List[Dict], metadata: Dict, output_file: Path,
                        offtarget_file: Path, specificity_file: Path,
                        min_size: int = 50, max_size: int = 10000):
    """Generate interactive HTML report"""

    # Did the user declare intent for any primer?
    references_provided = bool(metadata.get('references_provided'))
    references_per_primer = metadata.get('references_per_primer', {}) or {}

    # Summary statistics
    total_amplicons = len(results)
    targets = sum(1 for r in results if is_target(r))
    offtargets = sum(1 for r in results if is_off_target(r))
    unclassified = sum(1 for r in results if is_unclassified(r))

    # P(amplification) bins — always meaningful, computed from thermodynamics
    # regardless of intent declaration.
    high_risk = sum(1 for r in results if r['risk_classification'] == 'high')
    medium_risk = sum(1 for r in results if r['risk_classification'] == 'medium')
    low_risk = sum(1 for r in results if r['risk_classification'] == 'low')

    # Specificity / off-target sections only render when intent was declared.
    primer_specificity = []
    problem_primers = []
    if references_provided and specificity_file.exists():
        with open(specificity_file, 'r') as f:
            headers = f.readline().strip().split('\t')
            for line in f:
                if line.strip():
                    values = line.strip().split('\t')
                    primer_specificity.append(dict(zip(headers, values)))
        problem_primers = [p for p in primer_specificity
                          if int(p.get('High_Risk_OffTargets', 0)) > 0]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Thermodynamic Analysis Report</title>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            margin-bottom: 10px;
        }}
        .subtitle {{
            color: #666;
            margin-bottom: 30px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 6px;
            border-left: 4px solid #007bff;
        }}
        .stat-card.high {{ border-left-color: #dc3545; }}
        .stat-card.medium {{ border-left-color: #ffc107; }}
        .stat-card.low {{ border-left-color: #28a745; }}
        .stat-value {{
            font-size: 32px;
            font-weight: bold;
            color: #333;
        }}
        .stat-label {{
            color: #666;
            font-size: 14px;
            margin-top: 5px;
        }}
        .section {{
            margin: 40px 0;
        }}
        .section h2 {{
            color: #333;
            border-bottom: 2px solid #007bff;
            padding-bottom: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .risk-high {{ color: #dc3545; font-weight: bold; }}
        .risk-medium {{ color: #ffc107; font-weight: bold; }}
        .risk-low {{ color: #28a745; }}
        .metadata {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            margin: 20px 0;
            font-size: 14px;
        }}
        .metadata strong {{ color: #333; }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
            font-size: 12px;
            text-align: center;
        }}
        .rating-excellent {{ color: #28a745; font-weight: bold; }}
        .rating-good {{ color: #17a2b8; font-weight: bold; }}
        .rating-fair {{ color: #ffc107; font-weight: bold; }}
        .rating-poor {{ color: #dc3545; font-weight: bold; }}
        .dataTables_wrapper {{
            padding: 20px 0;
        }}
        .alignment-row {{
            display: none;
            background: #f8f9fa;
        }}
        .alignment-row.show {{
            display: table-row;
        }}
        .alignment-container {{
            padding: 15px;
            font-family: 'Courier New', monospace;
            font-size: 13px;
            line-height: 1.6;
            background: white;
            border-radius: 4px;
            margin: 10px;
        }}
        .alignment-header {{
            font-weight: bold;
            color: #333;
            margin-bottom: 10px;
        }}
        .alignment-text {{
            white-space: pre;
            color: #333;
        }}
        .match-good {{
            color: #28a745;
        }}
        .match-warn {{
            color: #ffc107;
        }}
        .match-bad {{
            color: #dc3545;
        }}
        .mismatch-3prime {{
            font-weight: bold;
            color: #dc3545;
            background: #ffe5e5;
        }}
        .expand-btn {{
            cursor: pointer;
            color: #007bff;
            text-decoration: underline;
            font-size: 12px;
        }}
        .expand-btn:hover {{
            color: #0056b3;
        }}
        .mismatch-indicator {{
            display: inline-block;
            width: 100%;
            height: 4px;
            background: #e9ecef;
            position: relative;
            margin-top: 5px;
        }}
        .mismatch-bar {{
            position: absolute;
            height: 100%;
            background: linear-gradient(to right, #ffc107, #dc3545);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Thermodynamic Primer Analysis Report</h1>
        <div class="subtitle">Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>

        {_intent_banner_html(references_provided, references_per_primer)}

        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{total_amplicons}</div>
                <div class="stat-label">Total Amplicons</div>
            </div>
            <div class="stat-card high">
                <div class="stat-value">{high_risk}</div>
                <div class="stat-label">High P(amp)</div>
            </div>
            <div class="stat-card medium">
                <div class="stat-value">{medium_risk}</div>
                <div class="stat-label">Medium P(amp)</div>
            </div>
            <div class="stat-card low">
                <div class="stat-value">{low_risk}</div>
                <div class="stat-label">Low P(amp)</div>
            </div>
            {"<div class='stat-card'><div class='stat-value'>" + str(targets) + "</div><div class='stat-label'>Intended Targets</div></div>" if references_provided else ""}
            {"<div class='stat-card'><div class='stat-value'>" + str(offtargets) + "</div><div class='stat-label'>Off-Targets</div></div>" if references_provided else ""}
        </div>

        <div class="metadata">
            <strong>PCR Conditions:</strong> {metadata.get('master_mix', 'Unknown')}<br>
            <strong>Annealing Temperature:</strong> {metadata.get('annealing_temp', 'N/A')}°C<br>
            <strong>Buffer Composition:</strong>
            Na⁺ = {metadata['pcr_conditions'].get('Na', 'N/A')} mM,
            Mg²⁺ = {metadata['pcr_conditions'].get('Mg', 'N/A')} mM,
            dNTPs = {metadata['pcr_conditions'].get('dNTPs', 'N/A')} mM<br>
            <strong>ΔTm Thresholds:</strong>
            High risk < {metadata['delta_tm_thresholds']['high_risk']}°C,
            Medium risk < {metadata['delta_tm_thresholds']['medium_risk']}°C<br>
            <strong>Amplicon Size Filter:</strong> {min_size:,} - {max_size:,} bp
        </div>

        <div class="section" style="background: #e3f2fd; border-left: 4px solid #2196F3; padding: 15px; margin: 20px 0;">
            <h3 style="margin-top: 0; color: #1565C0;">Understanding Alignment Results</h3>
            <p><strong>Degenerate Primers:</strong> Primers containing IUPAC ambiguity codes (e.g., [TCY] = T or C or Y)
            represent multiple variants. Primersearch matches any variant, but alignments show one specific variant.</p>

            <p><strong>Two Mismatch Counts:</strong></p>
            <ul style="margin: 5px 0;">
                <li><strong>Variant mismatches:</strong> Mismatches between the selected variant and genome sequence (shown in alignment)</li>
                <li><strong>Degenerate mismatches:</strong> Mismatches considering all possible primer variants (primersearch count)</li>
            </ul>

            <p><strong>Matched Variant:</strong> For each degenerate position, the variant that best matches the genome is selected.
            For intended targets with 0 degenerate mismatches, the alignment should show perfect matches.</p>

            <p><strong>3' End Mismatches:</strong> Mismatches within the last 5 bases from the 3' end are critical for PCR extension.
            Even 1 mismatch at the 3' end can significantly reduce amplification efficiency.</p>

            <p><strong>Gaps in Alignment:</strong> Multiple gaps (>3) suggest position extraction may be imprecise.
            This can occur at sequence boundaries or with complex primers.</p>
        </div>
"""

    # Add Primer Specificity section (only when intent was declared)
    if references_provided and primer_specificity:
        html += """
        <div class="section">
            <h2>Primer Specificity Summary</h2>
            <p>Overall performance metrics for each primer pair.</p>
            <table id="specificityTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Total Hits</th>
                        <th>Target Hits</th>
                        <th>High Risk Off-Targets</th>
                        <th>Specificity Score</th>
                        <th>Rating</th>
                    </tr>
                </thead>
                <tbody>
"""
        for p in primer_specificity:
            rating = p.get('Specificity_Rating', 'N/A')
            rating_class = f"rating-{rating.lower()}"
            html += f"""
                    <tr>
                        <td>{p.get('Primer', 'N/A')}</td>
                        <td>{p.get('Total_Hits', 'N/A')}</td>
                        <td>{p.get('Target_Hits', 'N/A')}</td>
                        <td>{p.get('High_Risk_OffTargets', 'N/A')}</td>
                        <td>{p.get('Specificity_Score', 'N/A')}</td>
                        <td class="{rating_class}">{rating}</td>
                    </tr>
"""
        html += """
                </tbody>
            </table>
        </div>
"""

    # Add Primer QC section (GC content and secondary structure)
    # Aggregate primer-level metrics (one entry per primer)
    primer_qc_data = {}
    for r in results:
        pname = r['primer_name']
        if pname not in primer_qc_data:
            primer_qc_data[pname] = {
                'fwd_gc': r.get('forward_gc_content', 'N/A'),
                'rev_gc': r.get('reverse_gc_content', 'N/A'),
                'fwd_dimer': r.get('forward_self_dimer_risk', 'N/A'),
                'rev_dimer': r.get('reverse_self_dimer_risk', 'N/A'),
                'fwd_hairpin': r.get('forward_hairpin_risk', 'N/A'),
                'rev_hairpin': r.get('reverse_hairpin_risk', 'N/A')
            }

    if primer_qc_data:
        html += """
        <div class="section">
            <h2>Primer Quality Control</h2>
            <p>GC content and secondary structure analysis for primer pairs.</p>
            <table id="primerQCTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Forward GC%</th>
                        <th>Reverse GC%</th>
                        <th>Fwd Dimer Risk</th>
                        <th>Rev Dimer Risk</th>
                        <th>Fwd Hairpin Risk</th>
                        <th>Rev Hairpin Risk</th>
                    </tr>
                </thead>
                <tbody>
"""
        for pname, qc in sorted(primer_qc_data.items()):
            # Color code risk levels
            fwd_dimer_class = 'risk-high' if qc['fwd_dimer'] == 'high' else ('risk-medium' if qc['fwd_dimer'] == 'medium' else 'risk-low')
            rev_dimer_class = 'risk-high' if qc['rev_dimer'] == 'high' else ('risk-medium' if qc['rev_dimer'] == 'medium' else 'risk-low')
            fwd_hairpin_class = 'risk-high' if qc['fwd_hairpin'] == 'high' else ('risk-medium' if qc['fwd_hairpin'] == 'medium' else 'risk-low')
            rev_hairpin_class = 'risk-high' if qc['rev_hairpin'] == 'high' else ('risk-medium' if qc['rev_hairpin'] == 'medium' else 'risk-low')

            html += f"""
                    <tr>
                        <td>{pname}</td>
                        <td>{qc['fwd_gc']}</td>
                        <td>{qc['rev_gc']}</td>
                        <td class="{fwd_dimer_class}">{qc['fwd_dimer']}</td>
                        <td class="{rev_dimer_class}">{qc['rev_dimer']}</td>
                        <td class="{fwd_hairpin_class}">{qc['fwd_hairpin']}</td>
                        <td class="{rev_hairpin_class}">{qc['rev_hairpin']}</td>
                    </tr>
"""
        html += """
                </tbody>
            </table>
        </div>
"""

    # Add Intended Target Validation section (only when intent declared)
    target_amplicons = [r for r in results if is_target(r)] if references_provided else []
    if target_amplicons:
        html += f"""
        <div class="section">
            <h2>Intended Target Validation</h2>
            <p>Verification that primers bind to their intended reference sequences ({len(target_amplicons)} targets detected).</p>
            <table id="targetsTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Position</th>
                        <th>Amplicon Length</th>
                        <th>Mismatches</th>
                        <th>ΔTm (°C)</th>
                    </tr>
                </thead>
                <tbody>
"""
        for r in target_amplicons:  # Show all
            # Use primersearch mismatch count for intended targets
            # (primersearch is IUPAC-aware, so this reflects what the user specified)
            mm_count = r.get('total_mismatches_primersearch', r['total_mismatches'])
            pos_str = format_position_range(
                r.get('forward_position'),
                r.get('reverse_position'),
                r.get('amplicon_length')
            )

            html += f"""
                    <tr>
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{pos_str}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{mm_count}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                    </tr>
"""
        html += """
                </tbody>
            </table>
        </div>
"""

    # Add Problem Primers section (only when intent declared)
    if references_provided and problem_primers:
        html += f"""
        <div class="section">
            <h2>Problem Primers</h2>
            <p>Primers with high off-target risk ({len(problem_primers)} primers flagged).</p>
            <table id="problemPrimersTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>High Risk Off-Targets</th>
                        <th>Medium Risk Off-Targets</th>
                        <th>Specificity Score</th>
                        <th>Recommendation</th>
                    </tr>
                </thead>
                <tbody>
"""
        for p in problem_primers:
            high_risk_count = int(p.get('High_Risk_OffTargets', 0))
            score = float(p.get('Specificity_Score', 0))
            if high_risk_count > 5:
                recommendation = "Consider redesign - Multiple high-risk off-targets"
            elif score < 75:
                recommendation = "Review primers - Low specificity"
            else:
                recommendation = "Monitor - Some off-target risk"

            html += f"""
                    <tr>
                        <td>{p.get('Primer', 'N/A')}</td>
                        <td class="risk-high">{p.get('High_Risk_OffTargets', 'N/A')}</td>
                        <td class="risk-medium">{p.get('Medium_Risk_OffTargets', 'N/A')}</td>
                        <td>{p.get('Specificity_Score', 'N/A')}</td>
                        <td>{recommendation}</td>
                    </tr>
"""
        html += """
                </tbody>
            </table>
        </div>
"""

    # Off-target sections are only rendered when intent was declared. Without
    # a reference for the primer, we can't say a hit is "off-target" — it's
    # simply unclassified, and forcing it into a risk table would be misleading.
    if references_provided:
        html += f"""
        <div class="section">
            <h2>High-Risk Off-Target Amplifications</h2>
            <p>Amplicons with high probability of unintended amplification (ΔTm < {metadata['delta_tm_thresholds']['high_risk']}°C and minimal 3' mismatches).</p>
"""

        high_risk_offtargets = [r for r in results
                               if is_off_target(r)
                               and r['risk_classification'] == 'high']

        if high_risk_offtargets:
            html += """
            <table id="highRiskTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Position</th>
                        <th>Length (bp)</th>
                        <th>P(Amplification)</th>
                        <th>ΔTm (°C)</th>
                        <th>Mismatches</th>
                        <th>Alignment</th>
                    </tr>
                </thead>
                <tbody>
"""
            for idx, r in enumerate(high_risk_offtargets):
                row_id = f"high_risk_{idx}"
                has_alignment = 'forward_alignment' in r or 'reverse_alignment' in r

                if has_alignment:
                    alignment_html = format_alignment_html(r, row_id)
                    alignment_data = alignment_html.replace('"', '&quot;').replace("'", '&#39;')
                    alignment_btn = f'<span class="expand-btn" onclick="toggleAlignment(\'{row_id}\', this)" data-alignment="{alignment_data}">Show</span>'
                else:
                    alignment_btn = 'N/A'

                pos_str = format_position_range(
                    r.get('forward_position'),
                    r.get('reverse_position'),
                    r.get('amplicon_length')
                )

                html += f"""
                    <tr data-row-id="{row_id}">
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{pos_str}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['p_amplification_overall']:.3f}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                        <td>{r['total_mismatches']}</td>
                        <td>{alignment_btn}</td>
                    </tr>
"""
            html += "                </tbody>\n            </table>\n"
        else:
            html += "            <p>No high-risk off-target amplifications detected.</p>\n"

        # Medium risk off-targets table
        html += """
        </div>

        <div class="section">
            <h2>Medium Risk Off-Targets</h2>
            <p>Off-target amplifications with moderate probability (0.3-0.7). These may amplify under some conditions.</p>
"""

        medium_risk_offtargets = [r for r in results
                                  if is_off_target(r)
                                  and r['risk_classification'] == 'medium']

        if medium_risk_offtargets:
            html += """
            <table id="mediumRiskTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Position</th>
                        <th>Length (bp)</th>
                        <th>P(Amplification)</th>
                        <th>ΔTm (°C)</th>
                        <th>Mismatches</th>
                        <th>Alignment</th>
                    </tr>
                </thead>
                <tbody>
"""
            for idx, r in enumerate(medium_risk_offtargets):
                row_id = f"medium_risk_{idx}"
                has_alignment = 'forward_alignment' in r or 'reverse_alignment' in r

                if has_alignment:
                    alignment_html = format_alignment_html(r, row_id)
                    alignment_data = alignment_html.replace('"', '&quot;').replace("'", '&#39;')
                    alignment_btn = f'<span class="expand-btn" onclick="toggleAlignment(\'{row_id}\', this)" data-alignment="{alignment_data}">Show</span>'
                else:
                    alignment_btn = 'N/A'

                pos_str = format_position_range(
                    r.get('forward_position'),
                    r.get('reverse_position'),
                    r.get('amplicon_length')
                )

                html += f"""
                    <tr data-row-id="{row_id}">
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{pos_str}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['p_amplification_overall']:.3f}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                        <td>{r['total_mismatches']}</td>
                        <td>{alignment_btn}</td>
                    </tr>
"""
            html += "                </tbody>\n            </table>\n"
        else:
            html += "            <p>No medium-risk off-target amplifications detected.</p>\n"

        # Low risk off-targets table
        html += """
        </div>

        <div class="section">
            <h2>Low Risk Off-Targets</h2>
            <p>Off-target amplifications with low probability (<0.3). These are unlikely to amplify significantly.</p>
"""

        low_risk_offtargets = [r for r in results
                               if is_off_target(r)
                               and r['risk_classification'] == 'low']

        if low_risk_offtargets:
            html += """
            <table id="lowRiskTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Position</th>
                        <th>Length (bp)</th>
                        <th>P(Amplification)</th>
                        <th>ΔTm (°C)</th>
                        <th>Mismatches</th>
                        <th>Alignment</th>
                    </tr>
                </thead>
                <tbody>
"""
            for idx, r in enumerate(low_risk_offtargets):
                row_id = f"low_risk_{idx}"
                has_alignment = 'forward_alignment' in r or 'reverse_alignment' in r

                if has_alignment:
                    alignment_html = format_alignment_html(r, row_id)
                    alignment_data = alignment_html.replace('"', '&quot;').replace("'", '&#39;')
                    alignment_btn = f'<span class="expand-btn" onclick="toggleAlignment(\'{row_id}\', this)" data-alignment="{alignment_data}">Show</span>'
                else:
                    alignment_btn = 'N/A'

                pos_str = format_position_range(
                    r.get('forward_position'),
                    r.get('reverse_position'),
                    r.get('amplicon_length')
                )

                html += f"""
                    <tr data-row-id="{row_id}">
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{pos_str}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['p_amplification_overall']:.3f}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                        <td>{r['total_mismatches']}</td>
                        <td>{alignment_btn}</td>
                    </tr>
"""
            html += "                </tbody>\n            </table>\n"
        else:
            html += "            <p>No low-risk off-target amplifications detected.</p>\n"

        html += "        </div>\n"

    html += f"""
        <div class="section">
            <h2>Files Generated</h2>
            <ul>
                <li><strong>thermo_analysis.tsv</strong> - Detailed per-amplicon thermodynamic metrics</li>
                <li><strong>offtarget_summary.tsv</strong> - Summary of medium/high-risk off-targets</li>
                <li><strong>primer_specificity.tsv</strong> - Per-primer specificity scores and statistics</li>
            </ul>
        </div>

        <div class="footer">
            Generated by Genome Primer Search Pipeline - Thermodynamic Analysis Module
        </div>
    </div>

    <script>
        $(document).ready(function() {{
            // Initialize DataTables with pagination and sorting
            $('#specificityTable').DataTable({{
                pageLength: 10,
                order: [[4, 'desc']]  // Sort by Specificity Score descending
            }});

            $('#primerQCTable').DataTable({{
                pageLength: 25,
                order: [[0, 'asc']]  // Sort by Primer name ascending
            }});

            $('#targetsTable').DataTable({{
                pageLength: 25,
                order: [[0, 'asc']]  // Sort by Primer name ascending
            }});

            $('#problemPrimersTable').DataTable({{
                pageLength: 10,
                order: [[1, 'desc']]  // Sort by High Risk Off-Targets descending
            }});

            $('#highRiskTable').DataTable({{
                pageLength: 25,
                order: [[3, 'desc']]  // Sort by P(Amplification) descending
            }});

            $('#mediumRiskTable').DataTable({{
                pageLength: 25,
                order: [[3, 'desc']]  // Sort by P(Amplification) descending
            }});

            $('#lowRiskTable').DataTable({{
                pageLength: 25,
                order: [[3, 'desc']]  // Sort by P(Amplification) descending
            }});
        }});

        // Toggle alignment visibility with dynamic row injection
        function toggleAlignment(rowId, btnElement) {{
            const parentRow = btnElement.closest('tr');
            const existingAlignmentRow = parentRow.nextElementSibling;

            // Check if alignment row already exists
            if (existingAlignmentRow && existingAlignmentRow.classList.contains('alignment-row')) {{
                // Remove it
                existingAlignmentRow.remove();
                btnElement.textContent = 'Show';
            }} else {{
                // Create and insert alignment row
                const alignmentData = btnElement.getAttribute('data-alignment');
                if (alignmentData) {{
                    const newRow = document.createElement('tr');
                    newRow.className = 'alignment-row show';
                    newRow.id = rowId;
                    newRow.innerHTML = '<td colspan="8">' + alignmentData + '</td>';

                    // Insert after current row
                    parentRow.parentNode.insertBefore(newRow, parentRow.nextSibling);
                    btnElement.textContent = 'Hide';
                }}
            }}
        }}
    </script>
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description='Generate thermodynamic analysis reports from JSON'
    )

    parser.add_argument('json_files', type=Path, nargs='+',
                       help='Input JSON file(s) from thermodynamic_analysis.py')
    parser.add_argument('-o', '--outdir', type=Path, required=True,
                       help='Output directory for reports')
    parser.add_argument('--min-size', type=int, default=50,
                       help='Minimum amplicon size filter (bp)')
    parser.add_argument('--max-size', type=int, default=10000,
                       help='Maximum amplicon size filter (bp)')

    args = parser.parse_args()

    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Aggregate results from all JSON files
    all_results = []
    metadata = None

    for json_file in args.json_files:
        print(f"Loading {json_file}", file=sys.stderr)
        with open(json_file) as f:
            data = json.load(f)

        if metadata is None:
            metadata = data['metadata']

        all_results.extend(data['results'])

    print(f"Loaded {len(all_results)} amplicons from {len(args.json_files)} file(s)", file=sys.stderr)

    # Generate reports
    detailed_file = args.outdir / 'thermo_analysis.tsv'
    offtarget_file = args.outdir / 'offtarget_summary.tsv'
    specificity_file = args.outdir / 'primer_specificity.tsv'
    html_file = args.outdir / 'thermo_analysis.html'

    print(f"Generating detailed TSV: {detailed_file}", file=sys.stderr)
    generate_detailed_tsv(all_results, detailed_file)

    print(f"Generating off-target summary: {offtarget_file}", file=sys.stderr)
    generate_offtarget_summary(all_results, offtarget_file)

    print(f"Generating primer specificity report: {specificity_file}", file=sys.stderr)
    generate_primer_specificity(all_results, specificity_file)

    print(f"Generating HTML report: {html_file}", file=sys.stderr)
    generate_html_report(all_results, metadata, html_file, offtarget_file, specificity_file,
                        args.min_size, args.max_size)

    print("\nReport generation complete!", file=sys.stderr)
    print(f"  Detailed analysis: {detailed_file}", file=sys.stderr)
    print(f"  Off-target summary: {offtarget_file}", file=sys.stderr)
    print(f"  Primer specificity: {specificity_file}", file=sys.stderr)
    print(f"  HTML report: {html_file}", file=sys.stderr)


if __name__ == '__main__':
    main()
