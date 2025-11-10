#!/usr/bin/env python3
"""
Filter amplicons from EMBOSS primersearch output based on size criteria.
Removes likely circular genome artifacts and invalid amplicons.
"""

import argparse
import json
import re
import sys


def parse_primersearch_output(input_file):
    """
    Parse EMBOSS primersearch output and extract amplicon information.
    Returns list of amplicons with metadata.
    """
    amplicons = []
    current_primer = None
    current_amplicon = {}
    in_amplicon = False

    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            # Detect primer name
            if line.startswith('Primer name '):
                current_primer = line.replace('Primer name ', '').strip()
                in_amplicon = False

            # Detect amplicon start
            elif line.startswith('Amplimer '):
                if current_amplicon and 'length' in current_amplicon:
                    amplicons.append(current_amplicon)

                current_amplicon = {
                    'primer': current_primer,
                    'header': line,
                    'lines': [line]
                }
                in_amplicon = True

            # Detect amplicon length
            elif in_amplicon and 'Amplimer length:' in line:
                match = re.search(r'Amplimer length:\s+(\d+)\s+bp', line)
                if match:
                    current_amplicon['length'] = int(match.group(1))
                    current_amplicon['lines'].append(line)

            # Collect all lines for current amplicon
            elif in_amplicon:
                current_amplicon['lines'].append(line)

    # Add last amplicon if exists
    if current_amplicon and 'length' in current_amplicon:
        amplicons.append(current_amplicon)

    return amplicons


def filter_amplicons(amplicons, min_size, max_size):
    """
    Filter amplicons based on size criteria.
    Returns filtered amplicons and statistics.
    """
    filtered = []
    stats = {
        'total_amplicons': len(amplicons),
        'retained': 0,
        'filtered_too_small': 0,
        'filtered_too_large': 0,
        'primers_with_hits': set(),
        'primers_after_filter': set()
    }

    for amp in amplicons:
        stats['primers_with_hits'].add(amp['primer'])
        length = amp['length']

        if length < min_size:
            stats['filtered_too_small'] += 1
        elif length > max_size:
            stats['filtered_too_large'] += 1
        else:
            filtered.append(amp)
            stats['retained'] += 1
            stats['primers_after_filter'].add(amp['primer'])

    # Convert sets to counts for JSON serialization
    stats['unique_primers_before'] = len(stats['primers_with_hits'])
    stats['unique_primers_after'] = len(stats['primers_after_filter'])
    stats['primers_with_hits'] = list(stats['primers_with_hits'])
    stats['primers_after_filter'] = list(stats['primers_after_filter'])

    return filtered, stats


def write_filtered_output(filtered_amplicons, output_file):
    """
    Write filtered amplicons back to primersearch format.
    """
    current_primer = None

    with open(output_file, 'w') as f:
        for amp in filtered_amplicons:
            # Write primer name header if changed
            if amp['primer'] != current_primer:
                if current_primer is not None:
                    f.write('\n')  # Blank line between primers
                f.write(f"Primer name {amp['primer']}\n")
                current_primer = amp['primer']

            # Write all amplicon lines
            for line in amp['lines']:
                f.write(line + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Filter amplicons from primersearch output'
    )
    parser.add_argument('--input', required=True, help='Input primersearch file')
    parser.add_argument('--output', required=True, help='Output filtered file')
    parser.add_argument('--stats', required=True, help='Output statistics JSON file')
    parser.add_argument('--min-size', type=int, default=50, help='Minimum amplicon size (bp)')
    parser.add_argument('--max-size', type=int, default=10000, help='Maximum amplicon size (bp)')

    args = parser.parse_args()

    # Parse input
    amplicons = parse_primersearch_output(args.input)

    # Filter amplicons
    filtered, stats = filter_amplicons(amplicons, args.min_size, args.max_size)

    # Write outputs
    write_filtered_output(filtered, args.output)

    with open(args.stats, 'w') as f:
        json.dump(stats, f, indent=2)

    # Print summary to stderr
    print(f"Filtered {stats['total_amplicons']} amplicons:", file=sys.stderr)
    print(f"  - Retained: {stats['retained']}", file=sys.stderr)
    print(f"  - Too small (<{args.min_size} bp): {stats['filtered_too_small']}", file=sys.stderr)
    print(f"  - Too large (>{args.max_size} bp): {stats['filtered_too_large']}", file=sys.stderr)
    print(f"  - Primers: {stats['unique_primers_before']} -> {stats['unique_primers_after']}", file=sys.stderr)


if __name__ == '__main__':
    main()
