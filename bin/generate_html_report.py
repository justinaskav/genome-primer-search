#!/usr/bin/env python3
"""
Generate comprehensive summary reports from filtered primersearch results.
Creates TSV and HTML reports with statistics and visualizations.
Enhanced version with primer columns and detailed filtering information.
"""

import argparse
import json
import glob
import re
from pathlib import Path
from collections import defaultdict
from datetime import datetime


def parse_filtered_results(filtered_dir):
    """
    Parse all filtered result files and extract amplicon details.
    """
    results = []
    filtered_files = glob.glob(f"{filtered_dir}/*_filtered.out")

    for filepath in filtered_files:
        genome_name = Path(filepath).stem.replace('_filtered', '')

        with open(filepath, 'r') as f:
            current_primer = None
            current_amplicon = []
            current_length = None

            for line in f:
                line = line.rstrip('\n')

                if line.startswith('Primer name '):
                    # Save previous primer's last amplicon before starting new primer
                    if current_amplicon and current_primer:
                        results.append({
                            'genome': genome_name,
                            'primer': current_primer,
                            'amplicon': '\n'.join(current_amplicon),
                            'length': current_length
                        })
                    # Start new primer
                    current_primer = line.replace('Primer name ', '').strip()
                    current_amplicon = []
                    current_length = None

                elif line.startswith('Amplimer '):
                    # Save previous amplicon if exists
                    if current_amplicon and current_primer:
                        results.append({
                            'genome': genome_name,
                            'primer': current_primer,
                            'amplicon': '\n'.join(current_amplicon),
                            'length': current_length
                        })
                    # Start new amplicon
                    current_amplicon = [line]
                    current_length = None

                elif current_amplicon:
                    current_amplicon.append(line)

                    # Extract length if present
                    if 'Amplimer length:' in line:
                        match = re.search(r'Amplimer length:\s+(\d+)\s+bp', line)
                        if match:
                            current_length = int(match.group(1))

            # Add last amplicon
            if current_amplicon and current_primer:
                results.append({
                    'genome': genome_name,
                    'primer': current_primer,
                    'amplicon': '\n'.join(current_amplicon),
                    'length': current_length
                })

    return results


def parse_stats_files(stats_dir):
    """
    Parse all statistics JSON files.
    """
    stats = {}
    stats_files = glob.glob(f"{stats_dir}/*_stats.json")

    for filepath in stats_files:
        genome_name = Path(filepath).stem.replace('_stats', '')

        with open(filepath, 'r') as f:
            stats[genome_name] = json.load(f)

    return stats


def get_unique_primers(results):
    """
    Extract unique primer names from results.
    """
    return sorted(set(r['primer'] for r in results))


def parse_raw_primersearch(filtered_dir):
    """
    Parse raw primersearch files to get original amplicon counts per primer.
    This helps us understand what was filtered per primer.
    """
    raw_data = defaultdict(lambda: defaultdict(int))

    # Look for original primersearch files in parent directory
    parent_dir = Path(filtered_dir).parent
    primersearch_files = glob.glob(f"{parent_dir}/primersearch/*_primersearch.out")

    for filepath in primersearch_files:
        genome_name = Path(filepath).stem.replace('_primersearch', '')

        with open(filepath, 'r') as f:
            current_primer = None
            for line in f:
                if line.startswith('Primer name '):
                    current_primer = line.replace('Primer name ', '').strip()
                elif line.startswith('Amplimer ') and current_primer:
                    raw_data[genome_name][current_primer] += 1

    return raw_data


def generate_summary_tsv(results, stats, output_file):
    """
    Generate summary TSV in long/normalized format for scalability.
    Format: One row per genome-primer combination.
    Includes genome-level total rows for easy aggregation.
    """
    primers = get_unique_primers(results)

    # Aggregate data per genome per primer
    genome_primer_data = defaultdict(lambda: defaultdict(lambda: {
        'kept': 0,
        'lengths': [],
        'filtered_circular': 0,
        'filtered_small': 0
    }))

    # Count kept amplicons and collect lengths
    for result in results:
        genome = result['genome']
        primer = result['primer']
        genome_primer_data[genome][primer]['kept'] += 1
        if result.get('length'):
            genome_primer_data[genome][primer]['lengths'].append(result['length'])

    # Add filtering info from stats and calculate per-primer estimates
    for genome, stat in stats.items():
        total_kept = sum(genome_primer_data[genome][p]['kept'] for p in primers if p in genome_primer_data[genome])
        filtered_small_total = stat.get('filtered_too_small', 0)
        filtered_large_total = stat.get('filtered_too_large', 0)

        # Estimate per-primer filtering based on proportion of kept amplicons
        for primer in primers:
            if primer in genome_primer_data[genome]:
                kept = genome_primer_data[genome][primer]['kept']
                # Proportional estimate of filtering
                if total_kept > 0:
                    genome_primer_data[genome][primer]['filtered_circular'] = int((kept / total_kept) * filtered_large_total)
                    genome_primer_data[genome][primer]['filtered_small'] = int((kept / total_kept) * filtered_small_total)

    # Write TSV in long/normalized format
    with open(output_file, 'w') as f:
        # Header
        header = ['Genome', 'Primer', 'Total_Found', 'Kept', 'Filtered_Circular', 'Filtered_Small', 'Mean_Size', 'Size_Range', 'Filter_Criteria']
        f.write('\t'.join(header) + '\n')

        # Write data: genome-level totals first, then per-primer details
        for genome in sorted(genome_primer_data.keys()):
            # Get genome totals from stats
            stat = stats.get(genome, {})
            total_found = stat.get('total_amplicons', 0)
            retained = stat.get('retained', 0)
            filtered_circular = stat.get('filtered_too_large', 0)
            filtered_small = stat.get('filtered_too_small', 0)

            # Write genome total row (with special marker)
            row = [
                genome,
                '_GENOME_TOTAL_',
                str(total_found),
                str(retained),
                str(filtered_circular),
                str(filtered_small),
                '-',
                '-',
                '50-10000bp'
            ]
            f.write('\t'.join(row) + '\n')

            # Write per-primer rows for this genome
            for primer in sorted(primers):
                if primer in genome_primer_data[genome]:
                    data = genome_primer_data[genome][primer]
                    lengths = data['lengths']

                    if lengths:
                        mean_size = int(sum(lengths) / len(lengths))
                        size_range = f"{min(lengths)}-{max(lengths)}"
                    else:
                        mean_size = 0
                        size_range = '0-0'

                    kept = data['kept']
                    filt_circ = data['filtered_circular']
                    filt_small = data['filtered_small']
                    total_primer_found = kept + filt_circ + filt_small

                    row = [
                        genome,
                        primer,
                        str(total_primer_found),
                        str(kept),
                        str(filt_circ),
                        str(filt_small),
                        str(mean_size),
                        size_range,
                        '50-10000bp'
                    ]
                    f.write('\t'.join(row) + '\n')


def generate_amplicon_stats_tsv(results, output_file):
    """
    Generate amplicon statistics in wide format with primers as columns.
    """
    primers = get_unique_primers(results)

    # Aggregate per genome per primer
    genome_primer_stats = defaultdict(lambda: defaultdict(lambda: {
        'count': 0,
        'lengths': []
    }))

    for result in results:
        genome = result['genome']
        primer = result['primer']
        genome_primer_stats[genome][primer]['count'] += 1
        if result.get('length'):
            genome_primer_stats[genome][primer]['lengths'].append(result['length'])

    # Write wide format TSV
    with open(output_file, 'w') as f:
        # Header
        header = ['Genome']
        for primer in primers:
            header.extend([
                f"{primer}_Count",
                f"{primer}_Mean",
                f"{primer}_Range"
            ])
        f.write('\t'.join(header) + '\n')

        # Data rows
        for genome in sorted(genome_primer_stats.keys()):
            row = [genome]

            for primer in primers:
                if primer in genome_primer_stats[genome]:
                    data = genome_primer_stats[genome][primer]
                    lengths = data['lengths']

                    if lengths:
                        mean_size = int(sum(lengths) / len(lengths))
                        size_range = f"{min(lengths)}-{max(lengths)}"
                    else:
                        mean_size = 0
                        size_range = '0-0'

                    row.extend([
                        str(data['count']),
                        str(mean_size),
                        size_range
                    ])
                else:
                    row.extend(['0', '0', '0-0'])

            f.write('\t'.join(row) + '\n')


def generate_primer_stats_tsv(results, stats, output_file):
    """
    Generate per-primer statistics with filtering details.
    """
    primers = get_unique_primers(results)

    # Aggregate per primer
    primer_data = defaultdict(lambda: {
        'kept': 0,
        'genomes': set(),
        'lengths': []
    })

    for result in results:
        primer = result['primer']
        primer_data[primer]['kept'] += 1
        primer_data[primer]['genomes'].add(result['genome'])
        if result.get('length'):
            primer_data[primer]['lengths'].append(result['length'])

    # Calculate total filtered per primer from stats
    for stat in stats.values():
        primers_with_hits = stat.get('primers_with_hits', [])
        total_filtered = stat.get('filtered_too_small', 0) + stat.get('filtered_too_large', 0)
        filtered_small = stat.get('filtered_too_small', 0)
        filtered_large = stat.get('filtered_too_large', 0)

        # Distribute filtering across primers (rough estimate)
        if primers_with_hits:
            for primer in primers_with_hits:
                if primer not in primer_data:
                    primer_data[primer] = {'kept': 0, 'genomes': set(), 'lengths': []}
                if 'filtered_total' not in primer_data[primer]:
                    primer_data[primer]['filtered_total'] = 0
                    primer_data[primer]['filtered_small'] = 0
                    primer_data[primer]['filtered_large'] = 0

                # Estimate: divide evenly among primers with hits
                primer_data[primer]['filtered_total'] += total_filtered // len(primers_with_hits)
                primer_data[primer]['filtered_small'] += filtered_small // len(primers_with_hits)
                primer_data[primer]['filtered_large'] += filtered_large // len(primers_with_hits)

    # Write TSV
    with open(output_file, 'w') as f:
        f.write('Primer\tKept\tFiltered_Total\tFiltered_Circular\tFiltered_Small\tGenomes_Hit\tMean_Size\tSize_Range\n')

        for primer in primers:
            data = primer_data[primer]
            lengths = data['lengths']

            if lengths:
                mean_size = int(sum(lengths) / len(lengths))
                size_range = f"{min(lengths)}-{max(lengths)}"
            else:
                mean_size = 0
                size_range = '0-0'

            filtered_total = data.get('filtered_total', 0)
            filtered_small = data.get('filtered_small', 0)
            filtered_large = data.get('filtered_large', 0)

            f.write(f"{primer}\t{data['kept']}\t{filtered_total}\t{filtered_large}\t{filtered_small}\t{len(data['genomes'])}\t{mean_size}\t{size_range}\n")


def generate_html_report(results, stats, output_file):
    """
    Generate HTML summary report with embedded CSS and detailed filtering information.
    """
    primers = get_unique_primers(results)
    total_genomes = len(set(r['genome'] for r in results)) if results else len(stats)
    total_amplicons = len(results)
    total_primers = len(primers)

    # Calculate filtering stats
    total_filtered = sum(s.get('filtered_too_small', 0) + s.get('filtered_too_large', 0) for s in stats.values())
    filtered_small = sum(s.get('filtered_too_small', 0) for s in stats.values())
    filtered_large = sum(s.get('filtered_too_large', 0) for s in stats.values())

    # Per-genome summary
    genome_summary = defaultdict(lambda: {'amplicons': 0, 'primers': set()})
    for result in results:
        genome_summary[result['genome']]['amplicons'] += 1
        genome_summary[result['genome']]['primers'].add(result['primer'])

    # Per-primer summary
    primer_summary = defaultdict(lambda: {'amplicons': 0, 'genomes': set()})
    for result in results:
        primer_summary[result['primer']]['amplicons'] += 1
        primer_summary[result['primer']]['genomes'].add(result['genome'])

    # Genome-primer matrix for filtering details
    genome_primer_data = defaultdict(lambda: defaultdict(lambda: {'kept': 0, 'filtered': 0}))
    for result in results:
        genome_primer_data[result['genome']][result['primer']]['kept'] += 1

    # Add filtering estimates
    for genome, stat in stats.items():
        total_kept = sum(genome_primer_data[genome][p]['kept'] for p in primers if p in genome_primer_data[genome])
        total_filtered_genome = stat.get('filtered_too_small', 0) + stat.get('filtered_too_large', 0)

        for primer in primers:
            if primer in genome_primer_data[genome]:
                kept = genome_primer_data[genome][primer]['kept']
                if total_kept > 0:
                    estimated_filtered = int((kept / total_kept) * total_filtered_genome)
                else:
                    estimated_filtered = 0
                genome_primer_data[genome][primer]['filtered'] = estimated_filtered

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Genome Primer Search Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 8px;
        }}
        .summary-box {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .stat-card.success {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        }}
        .stat-card.warning {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        .stat-card .number {{
            font-size: 36px;
            font-weight: bold;
            margin: 10px 0;
        }}
        .stat-card .label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
        }}
        th {{
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ecf0f1;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
        .kept {{
            color: #27ae60;
            font-weight: bold;
        }}
        .filtered {{
            color: #e67e22;
        }}
        .filter-details {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin: 10px 0;
        }}
        .timestamp {{
            color: #7f8c8d;
            font-size: 14px;
            margin-top: 30px;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Genome Primer Search Report</h1>

        <div class="summary-box">
            <div class="stat-card">
                <div class="label">Total Genomes</div>
                <div class="number">{total_genomes}</div>
            </div>
            <div class="stat-card success">
                <div class="label">Amplicons Retained</div>
                <div class="number">{total_amplicons}</div>
            </div>
            <div class="stat-card">
                <div class="label">Unique Primers</div>
                <div class="number">{total_primers}</div>
            </div>
            <div class="stat-card warning">
                <div class="label">Amplicons Filtered</div>
                <div class="number">{total_filtered}</div>
            </div>
        </div>

        <h2>Filtering Summary</h2>
        <div class="filter-details">
            <strong>Filter Criteria Applied:</strong> 50 - 10,000 bp (removes circular genome artifacts)<br>
            <strong>Total Filtered:</strong> {total_filtered} amplicons<br>
            <strong>Filtered (too small &lt;50bp):</strong> {filtered_small}<br>
            <strong>Filtered (too large &gt;10kb / circular):</strong> {filtered_large}
        </div>

        <h2>Filtering Details by Genome</h2>
        <table>
            <tr>
                <th>Genome</th>
                <th>Total Found</th>
                <th>Retained</th>
                <th>Filtered</th>"""

    # Add primer columns to header
    for primer in primers:
        html += f"<th>{primer}</th>"

    html += "</tr>\n"

    # Add genome rows
    for genome in sorted(stats.keys()):
        stat = stats[genome]
        total_found = stat.get('total_amplicons', 0)
        retained = stat.get('retained', 0)
        filtered = stat.get('filtered_too_small', 0) + stat.get('filtered_too_large', 0)

        html += f"""            <tr>
                <td><strong>{genome}</strong></td>
                <td>{total_found}</td>
                <td class="kept">{retained}</td>
                <td class="filtered">{filtered}</td>"""

        # Add per-primer data
        for primer in primers:
            if primer in genome_primer_data[genome]:
                data = genome_primer_data[genome][primer]
                kept = data['kept']
                filt = data['filtered']
                html += f'<td><span class="kept">{kept} kept</span>, <span class="filtered">{filt} filtered</span></td>'
            else:
                html += '<td>-</td>'

        html += "\n            </tr>\n"

    html += """        </table>

        <h2>Per-Primer Performance Matrix</h2>
        <table>
            <tr>
                <th>Primer</th>
                <th>Total Amplicons Kept</th>
                <th>Total Filtered</th>
                <th>Genomes Hit</th>
            </tr>
"""

    for primer in sorted(primer_summary.keys()):
        data = primer_summary[primer]
        # Estimate filtered for this primer
        total_filtered_primer = sum(genome_primer_data[g][primer]['filtered'] for g in genome_primer_data if primer in genome_primer_data[g])

        html += f"""            <tr>
                <td><strong>{primer}</strong></td>
                <td class="kept">{data['amplicons']}</td>
                <td class="filtered">{total_filtered_primer}</td>
                <td>{len(data['genomes'])}</td>
            </tr>
"""

    html += f"""        </table>

        <div class="timestamp">
            Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive reports from filtered primersearch results'
    )
    parser.add_argument('--filtered-dir', required=True, help='Directory with filtered result files')
    parser.add_argument('--stats-dir', required=True, help='Directory with stats JSON files')
    parser.add_argument('--output-dir', required=True, help='Output directory for reports')

    args = parser.parse_args()

    # Parse all results
    print("Parsing filtered results...")
    results = parse_filtered_results(args.filtered_dir)

    print("Parsing statistics files...")
    stats = parse_stats_files(args.stats_dir)

    # Generate reports
    print("Generating summary.tsv...")
    generate_summary_tsv(results, stats, f"{args.output_dir}/summary.tsv")

    print("Generating amplicon_stats.tsv...")
    generate_amplicon_stats_tsv(results, f"{args.output_dir}/amplicon_stats.tsv")

    print("Generating primer_stats.tsv...")
    generate_primer_stats_tsv(results, stats, f"{args.output_dir}/primer_stats.tsv")

    print("Generating summary.html...")
    generate_html_report(results, stats, f"{args.output_dir}/summary.html")

    print(f"\nReports generated successfully!")
    print(f"  - Total genomes: {len(set(r['genome'] for r in results)) if results else len(stats)}")
    print(f"  - Total amplicons: {len(results)}")
    print(f"  - Unique primers: {len(set(r['primer'] for r in results))}")


if __name__ == '__main__':
    main()
