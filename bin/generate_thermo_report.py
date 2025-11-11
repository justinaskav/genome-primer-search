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


def generate_detailed_tsv(results: List[Dict], output_file: Path):
    """Generate detailed per-amplicon TSV report"""

    headers = [
        'Primer', 'Genome', 'Amplicon_Length',
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
                r['amplicon_length'],
                'Yes' if r['is_intended_target'] else 'No',
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
        'Primer', 'Genome', 'Amplicon_Length',
        'Risk_Classification', 'P_Amplification',
        'Delta_Tm_Average', 'Total_Mismatches',
        'Forward_3prime_MM', 'Reverse_3prime_MM',
        'Recommendation'
    ]

    # Filter to non-target hits with medium/high risk
    offtargets = [r for r in results
                  if not r['is_intended_target']
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
    """Generate per-primer specificity metrics"""

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

        if r['is_intended_target']:
            stats['target_hits'] += 1
        else:
            stats['offtarget_hits'] += 1
            stats['mean_delta_tm'].append(r['delta_tm_average'])
            stats['mean_p_amp_offtarget'].append(r['p_amplification_overall'])

            if r['risk_classification'] == 'high':
                stats['high_risk_offtargets'] += 1
            elif r['risk_classification'] == 'medium':
                stats['medium_risk_offtargets'] += 1
            else:
                stats['low_risk_offtargets'] += 1

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


def generate_html_report(results: List[Dict], metadata: Dict, output_file: Path,
                        offtarget_file: Path, specificity_file: Path):
    """Generate interactive HTML report"""

    # Summary statistics
    total_amplicons = len(results)
    targets = sum(1 for r in results if r['is_intended_target'])
    offtargets = total_amplicons - targets

    # Risk classifications (only for off-targets, not intended targets)
    high_risk = sum(1 for r in results
                   if not r['is_intended_target'] and r['risk_classification'] == 'high')
    medium_risk = sum(1 for r in results
                     if not r['is_intended_target'] and r['risk_classification'] == 'medium')
    low_risk = sum(1 for r in results
                  if not r['is_intended_target'] and r['risk_classification'] == 'low')

    # Read primer specificity data
    primer_specificity = []
    if specificity_file.exists():
        with open(specificity_file, 'r') as f:
            headers = f.readline().strip().split('\t')
            for line in f:
                if line.strip():
                    values = line.strip().split('\t')
                    primer_specificity.append(dict(zip(headers, values)))

    # Calculate ΔTm distribution for off-targets
    offtarget_delta_tms = [r['delta_tm_average'] for r in results
                           if not r['is_intended_target'] and r['delta_tm_average'] > 0]

    # Identify problem primers (high off-target rate)
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
    </style>
</head>
<body>
    <div class="container">
        <h1>Thermodynamic Primer Analysis Report</h1>
        <div class="subtitle">Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>

        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{total_amplicons}</div>
                <div class="stat-label">Total Amplicons</div>
            </div>
            <div class="stat-card high">
                <div class="stat-value">{high_risk}</div>
                <div class="stat-label">High Risk</div>
            </div>
            <div class="stat-card medium">
                <div class="stat-value">{medium_risk}</div>
                <div class="stat-label">Medium Risk</div>
            </div>
            <div class="stat-card low">
                <div class="stat-value">{low_risk}</div>
                <div class="stat-label">Low Risk</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{targets}</div>
                <div class="stat-label">Intended Targets</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{offtargets}</div>
                <div class="stat-label">Off-Targets</div>
            </div>
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
            Medium risk < {metadata['delta_tm_thresholds']['medium_risk']}°C
        </div>
"""

    # Add Primer Specificity section
    if primer_specificity:
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

    # Add Intended Target Validation section
    target_amplicons = [r for r in results if r['is_intended_target']]
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
                        <th>Amplicon Length</th>
                        <th>Mismatches</th>
                        <th>ΔTm (°C)</th>
                    </tr>
                </thead>
                <tbody>
"""
        for r in target_amplicons:  # Show all
            html += f"""
                    <tr>
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['total_mismatches']}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                    </tr>
"""
        html += """
                </tbody>
            </table>
        </div>
"""

    # Add Problem Primers section
    if problem_primers:
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

    html += f"""
        <div class="section">
            <h2>High-Risk Off-Target Amplifications</h2>
            <p>Amplicons with high probability of unintended amplification (ΔTm < {metadata['delta_tm_thresholds']['high_risk']}°C and minimal 3' mismatches).</p>
"""

    # Add high-risk off-targets table
    high_risk_offtargets = [r for r in results
                           if not r['is_intended_target']
                           and r['risk_classification'] == 'high']

    if high_risk_offtargets:
        html += """
            <table id="highRiskTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Length (bp)</th>
                        <th>P(Amplification)</th>
                        <th>ΔTm (°C)</th>
                        <th>Mismatches</th>
                    </tr>
                </thead>
                <tbody>
"""
        for r in high_risk_offtargets:  # Show all
            html += f"""
                    <tr>
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['p_amplification_overall']:.3f}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                        <td>{r['total_mismatches']}</td>
                    </tr>
"""
        html += "                </tbody>\n            </table>\n"
    else:
        html += "            <p>No high-risk off-target amplifications detected.</p>\n"

    # Add medium risk off-targets table
    html += """
        </div>

        <div class="section">
            <h2>Medium Risk Off-Targets</h2>
            <p>Off-target amplifications with moderate probability (0.3-0.7). These may amplify under some conditions.</p>
"""

    medium_risk_offtargets = [r for r in results
                              if not r['is_intended_target']
                              and r['risk_classification'] == 'medium']

    if medium_risk_offtargets:
        html += """
            <table id="mediumRiskTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Length (bp)</th>
                        <th>P(Amplification)</th>
                        <th>ΔTm (°C)</th>
                        <th>Mismatches</th>
                    </tr>
                </thead>
                <tbody>
"""
        for r in medium_risk_offtargets:
            html += f"""
                    <tr>
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['p_amplification_overall']:.3f}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                        <td>{r['total_mismatches']}</td>
                    </tr>
"""
        html += "                </tbody>\n            </table>\n"
    else:
        html += "            <p>No medium-risk off-target amplifications detected.</p>\n"

    # Add low risk off-targets table
    html += """
        </div>

        <div class="section">
            <h2>Low Risk Off-Targets</h2>
            <p>Off-target amplifications with low probability (<0.3). These are unlikely to amplify significantly.</p>
"""

    low_risk_offtargets = [r for r in results
                           if not r['is_intended_target']
                           and r['risk_classification'] == 'low']

    if low_risk_offtargets:
        html += """
            <table id="lowRiskTable" class="display">
                <thead>
                    <tr>
                        <th>Primer</th>
                        <th>Genome</th>
                        <th>Length (bp)</th>
                        <th>P(Amplification)</th>
                        <th>ΔTm (°C)</th>
                        <th>Mismatches</th>
                    </tr>
                </thead>
                <tbody>
"""
        for r in low_risk_offtargets:
            html += f"""
                    <tr>
                        <td>{r['primer_name']}</td>
                        <td>{r['genome']}</td>
                        <td>{r['amplicon_length']}</td>
                        <td>{r['p_amplification_overall']:.3f}</td>
                        <td>{r['delta_tm_average']:.1f}</td>
                        <td>{r['total_mismatches']}</td>
                    </tr>
"""
        html += "                </tbody>\n            </table>\n"
    else:
        html += "            <p>No low-risk off-target amplifications detected.</p>\n"

    html += f"""
        </div>

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
    generate_html_report(all_results, metadata, html_file, offtarget_file, specificity_file)

    print("\nReport generation complete!", file=sys.stderr)
    print(f"  Detailed analysis: {detailed_file}", file=sys.stderr)
    print(f"  Off-target summary: {offtarget_file}", file=sys.stderr)
    print(f"  Primer specificity: {specificity_file}", file=sys.stderr)
    print(f"  HTML report: {html_file}", file=sys.stderr)


if __name__ == '__main__':
    main()
