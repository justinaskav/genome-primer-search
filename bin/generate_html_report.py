#!/usr/bin/env python3
"""
Generate summary reports from filtered primersearch results.

Consumes the per-genome stats.json files produced by filter_amplicons.py
(which now contain the full kept-amplicon records, exact per-primer counts,
and per-contig breakdown — no more proportional allocation).
"""

import argparse
import glob
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path


def load_stats(stats_dir):
    """Return {genome_name: stats_dict} loaded from *_stats.json files."""
    out = {}
    for filepath in sorted(glob.glob(f"{stats_dir}/*_stats.json")):
        genome_name = Path(filepath).stem.replace('_stats', '')
        with open(filepath) as f:
            out[genome_name] = json.load(f)
    return out


def all_amplicons(stats):
    """Flatten kept amplicons across genomes, annotating each with its genome."""
    for genome, s in stats.items():
        for amp in s.get('amplicons', []):
            amp = dict(amp)
            amp['genome'] = genome
            yield amp


def unique_primers(stats):
    primers = set()
    for s in stats.values():
        primers.update(s.get('per_primer', {}).keys())
    return sorted(primers)


def format_position(amp):
    fwd = amp.get('forward') or {}
    rev = amp.get('reverse') or {}
    fp = fwd.get('pos')
    rp = rev.get('pos')
    if fp is not None and rp is not None:
        return f"{fp:,}..{rp:,}"
    if fp is not None:
        return f"{fp:,}.."
    if rp is not None:
        return f"..{rp:,}"
    return "N/A"


def format_mismatches(amp):
    fwd = (amp.get('forward') or {}).get('mismatches')
    rev = (amp.get('reverse') or {}).get('mismatches')
    if fwd is None and rev is None:
        return "N/A"
    return f"{fwd if fwd is not None else '?'}f / {rev if rev is not None else '?'}r"


def generate_summary_tsv(stats, output_file, min_size, max_size):
    """Long-format TSV: one row per (genome, primer), plus _GENOME_TOTAL_ rows."""
    primers = unique_primers(stats)

    with open(output_file, 'w') as f:
        f.write('\t'.join([
            'Genome', 'Primer', 'Total_Found', 'Kept',
            'Filtered_Too_Large', 'Filtered_Too_Small',
            'Mean_Size', 'Size_Range', 'Filter_Criteria'
        ]) + '\n')

        for genome in sorted(stats.keys()):
            s = stats[genome]
            totals = s.get('totals', {})

            f.write('\t'.join([
                genome, '_GENOME_TOTAL_',
                str(totals.get('total_amplicons', 0)),
                str(totals.get('retained', 0)),
                str(totals.get('filtered_too_large', 0)),
                str(totals.get('filtered_too_small', 0)),
                '-', '-', f'{min_size}-{max_size}bp'
            ]) + '\n')

            # Per-primer length stats from kept amplicons
            lengths_by_primer = defaultdict(list)
            for amp in s.get('amplicons', []):
                if amp.get('length') is not None:
                    lengths_by_primer[amp.get('primer')].append(amp['length'])

            for primer in primers:
                pp = s.get('per_primer', {}).get(primer)
                if not pp:
                    continue
                lengths = lengths_by_primer.get(primer, [])
                if lengths:
                    mean_size = int(sum(lengths) / len(lengths))
                    size_range = f"{min(lengths)}-{max(lengths)}"
                else:
                    mean_size = 0
                    size_range = '0-0'

                f.write('\t'.join([
                    genome, primer,
                    str(pp.get('found', 0)),
                    str(pp.get('kept', 0)),
                    str(pp.get('filtered_too_large', 0)),
                    str(pp.get('filtered_too_small', 0)),
                    str(mean_size), size_range,
                    f'{min_size}-{max_size}bp'
                ]) + '\n')


def generate_amplicon_stats_tsv(stats, output_file):
    """Wide format: one row per genome, primer-count columns."""
    primers = unique_primers(stats)

    with open(output_file, 'w') as f:
        header = ['Genome']
        for p in primers:
            header.extend([f"{p}_Count", f"{p}_Mean", f"{p}_Range"])
        f.write('\t'.join(header) + '\n')

        for genome in sorted(stats.keys()):
            s = stats[genome]
            row = [genome]
            lengths_by_primer = defaultdict(list)
            for amp in s.get('amplicons', []):
                if amp.get('length') is not None:
                    lengths_by_primer[amp.get('primer')].append(amp['length'])

            for p in primers:
                pp = s.get('per_primer', {}).get(p, {})
                count = pp.get('kept', 0)
                lengths = lengths_by_primer.get(p, [])
                if lengths:
                    mean_size = int(sum(lengths) / len(lengths))
                    size_range = f"{min(lengths)}-{max(lengths)}"
                else:
                    mean_size, size_range = 0, '0-0'
                row.extend([str(count), str(mean_size), size_range])
            f.write('\t'.join(row) + '\n')


def generate_primer_stats_tsv(stats, output_file):
    primers = unique_primers(stats)

    # Aggregate across genomes
    agg = {p: {'kept': 0, 'filtered_too_small': 0, 'filtered_too_large': 0,
               'genomes': set(), 'lengths': []} for p in primers}

    for genome, s in stats.items():
        for p, pp in s.get('per_primer', {}).items():
            if pp.get('kept', 0) > 0:
                agg[p]['genomes'].add(genome)
            agg[p]['kept'] += pp.get('kept', 0)
            agg[p]['filtered_too_small'] += pp.get('filtered_too_small', 0)
            agg[p]['filtered_too_large'] += pp.get('filtered_too_large', 0)
        for amp in s.get('amplicons', []):
            if amp.get('length') is not None:
                agg[amp.get('primer')]['lengths'].append(amp['length'])

    with open(output_file, 'w') as f:
        f.write('Primer\tKept\tFiltered_Total\tFiltered_Too_Large\tFiltered_Too_Small\tGenomes_Hit\tMean_Size\tSize_Range\n')
        for p in primers:
            d = agg[p]
            lengths = d['lengths']
            if lengths:
                mean_size = int(sum(lengths) / len(lengths))
                size_range = f"{min(lengths)}-{max(lengths)}"
            else:
                mean_size, size_range = 0, '0-0'
            filtered_total = d['filtered_too_large'] + d['filtered_too_small']
            f.write(f"{p}\t{d['kept']}\t{filtered_total}\t{d['filtered_too_large']}\t"
                    f"{d['filtered_too_small']}\t{len(d['genomes'])}\t{mean_size}\t{size_range}\n")


HTML_CSS = """
body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
.container { max-width: 1400px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
h2 { color: #34495e; margin-top: 30px; border-bottom: 2px solid #ecf0f1; padding-bottom: 8px; }
.summary-box { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }
.stat-card { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
.stat-card.success { background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); }
.stat-card.warning { background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); }
.stat-card .number { font-size: 36px; font-weight: bold; margin: 10px 0; }
.stat-card .label { font-size: 14px; opacity: 0.9; }
table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 14px; }
th { background-color: #3498db; color: white; padding: 12px; text-align: left; font-weight: 600; }
td { padding: 10px 12px; border-bottom: 1px solid #ecf0f1; }
tr:hover { background-color: #f8f9fa; }
.kept { color: #27ae60; font-weight: bold; }
.filtered { color: #e67e22; }
.mismatch { color: #c0392b; font-weight: bold; }
.zero-mm { color: #7f8c8d; }
.filter-details { background-color: #ecf0f1; padding: 15px; border-radius: 5px; margin: 10px 0; }
.timestamp { color: #7f8c8d; font-size: 14px; margin-top: 30px; text-align: center; }
"""


def generate_html_report(stats, output_file, min_size, max_size):
    primers = unique_primers(stats)
    total_genomes = len(stats)

    amplicons = list(all_amplicons(stats))
    total_amplicons = len(amplicons)
    total_primers = len(primers)

    totals_across = {
        'total': sum(s.get('totals', {}).get('total_amplicons', 0) for s in stats.values()),
        'kept': sum(s.get('totals', {}).get('retained', 0) for s in stats.values()),
        'small': sum(s.get('totals', {}).get('filtered_too_small', 0) for s in stats.values()),
        'large': sum(s.get('totals', {}).get('filtered_too_large', 0) for s in stats.values()),
    }
    totals_across['filtered'] = totals_across['small'] + totals_across['large']

    parts = []
    parts.append(f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Genome Primer Search Report</title>
<style>{HTML_CSS}</style></head><body><div class="container">
<h1>Genome Primer Search Report</h1>
<div class="summary-box">
  <div class="stat-card"><div class="label">Total Genomes</div><div class="number">{total_genomes}</div></div>
  <div class="stat-card success"><div class="label">Amplicons Retained</div><div class="number">{total_amplicons}</div></div>
  <div class="stat-card"><div class="label">Unique Primers</div><div class="number">{total_primers}</div></div>
  <div class="stat-card warning"><div class="label">Amplicons Filtered</div><div class="number">{totals_across['filtered']}</div></div>
</div>
<h2>Filtering Summary</h2>
<div class="filter-details">
  <strong>Size window:</strong> {min_size:,}–{max_size:,} bp<br>
  <strong>Total filtered:</strong> {totals_across['filtered']} amplicons<br>
  <strong>Filtered too small (&lt;{min_size:,} bp):</strong> {totals_across['small']}<br>
  <strong>Filtered too large (&gt;{max_size:,} bp):</strong> {totals_across['large']}<br>
  <em>Large amplicons typically represent cross-locus pairings (e.g. primers binding to multiple rrn operons treated as linear amplicons by primersearch) rather than true PCR products.</em>
</div>
""")

    parts.append("""<h2>Filtering Details by Genome</h2>
<table><tr>
<th>Genome</th><th>Total Found</th><th>Retained</th><th>Filtered (large)</th><th>Filtered (small)</th>""")
    for p in primers:
        parts.append(f"<th>{p}</th>")
    parts.append("</tr>\n")

    for genome in sorted(stats.keys()):
        s = stats[genome]
        t = s.get('totals', {})
        parts.append(f"""<tr>
<td><strong>{genome}</strong></td>
<td>{t.get('total_amplicons', 0)}</td>
<td class="kept">{t.get('retained', 0)}</td>
<td class="filtered">{t.get('filtered_too_large', 0)}</td>
<td class="filtered">{t.get('filtered_too_small', 0)}</td>""")
        for p in primers:
            pp = s.get('per_primer', {}).get(p)
            if pp:
                kept = pp.get('kept', 0)
                filt = pp.get('filtered_too_large', 0) + pp.get('filtered_too_small', 0)
                parts.append(f'<td><span class="kept">{kept} kept</span>, <span class="filtered">{filt} filtered</span></td>')
            else:
                parts.append('<td>-</td>')
        parts.append("</tr>\n")
    parts.append("</table>\n")

    # Per-primer aggregate
    parts.append("""<h2>Per-Primer Performance Matrix</h2>
<table><tr>
<th>Primer</th><th>Amplicons Kept</th><th>Total Filtered</th><th>Filtered (large)</th><th>Filtered (small)</th><th>Genomes Hit</th>
</tr>""")
    for p in primers:
        kept = 0
        f_large = 0
        f_small = 0
        genomes_hit = 0
        for genome, s in stats.items():
            pp = s.get('per_primer', {}).get(p)
            if not pp:
                continue
            kept += pp.get('kept', 0)
            f_large += pp.get('filtered_too_large', 0)
            f_small += pp.get('filtered_too_small', 0)
            if pp.get('kept', 0) > 0:
                genomes_hit += 1
        parts.append(f"""<tr>
<td><strong>{p}</strong></td>
<td class="kept">{kept}</td>
<td class="filtered">{f_large + f_small}</td>
<td class="filtered">{f_large}</td>
<td class="filtered">{f_small}</td>
<td>{genomes_hit}</td>
</tr>""")
    parts.append("</table>\n")

    # Detailed amplicons — now includes contig + mismatches
    parts.append("""<h2>Detailed Amplicon Results</h2>
<p>Kept amplicons with primer-binding positions, mismatch counts (forward / reverse), and source contig.</p>
<table><tr>
<th>Genome</th><th>Contig</th><th>Primer</th><th>Position (fwd..rev)</th><th>Length (bp)</th><th>Mismatches</th>
</tr>""")
    for amp in sorted(amplicons, key=lambda a: (a['genome'], a.get('sequence_id') or '', a.get('primer') or '', a.get('forward', {}).get('pos') or 0)):
        mm_str = format_mismatches(amp)
        fwd_mm = (amp.get('forward') or {}).get('mismatches', 0) or 0
        rev_mm = (amp.get('reverse') or {}).get('mismatches', 0) or 0
        mm_class = 'mismatch' if (fwd_mm + rev_mm) > 0 else 'zero-mm'
        parts.append(f"""<tr>
<td>{amp['genome']}</td>
<td>{amp.get('sequence_id') or '-'}</td>
<td>{amp.get('primer') or '-'}</td>
<td>{format_position(amp)}</td>
<td>{amp.get('length', '-')}</td>
<td class="{mm_class}">{mm_str}</td>
</tr>""")
    parts.append("</table>\n")

    parts.append(f"""<div class="timestamp">Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
</div></body></html>""")

    with open(output_file, 'w') as f:
        f.write(''.join(parts))


def main():
    parser = argparse.ArgumentParser(description='Generate reports from filtered primersearch results')
    parser.add_argument('--stats-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--min-size', type=int, default=50)
    parser.add_argument('--max-size', type=int, default=10000)
    # Accept --filtered-dir for backward compat with the existing Nextflow process; unused.
    parser.add_argument('--filtered-dir', required=False, default=None)
    args = parser.parse_args()

    print("Loading stats...")
    stats = load_stats(args.stats_dir)

    print("Generating summary.tsv...")
    generate_summary_tsv(stats, f"{args.output_dir}/summary.tsv", args.min_size, args.max_size)

    print("Generating amplicon_stats.tsv...")
    generate_amplicon_stats_tsv(stats, f"{args.output_dir}/amplicon_stats.tsv")

    print("Generating primer_stats.tsv...")
    generate_primer_stats_tsv(stats, f"{args.output_dir}/primer_stats.tsv")

    print("Generating summary.html...")
    generate_html_report(stats, f"{args.output_dir}/summary.html", args.min_size, args.max_size)

    n_amp = sum(len(s.get('amplicons', [])) for s in stats.values())
    n_primers = len(unique_primers(stats))
    print(f"\nReports generated successfully!")
    print(f"  - Genomes: {len(stats)}")
    print(f"  - Kept amplicons: {n_amp}")
    print(f"  - Unique primers: {n_primers}")


if __name__ == '__main__':
    main()
