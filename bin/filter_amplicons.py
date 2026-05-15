#!/usr/bin/env python3
"""
Filter amplicons from EMBOSS primersearch output based on size criteria.

Output schema (stats.json):
    {
      "min_size": 50,
      "max_size": 10000,
      "totals": {
        "total_amplicons": int,
        "retained": int,
        "filtered_too_small": int,
        "filtered_too_large": int
      },
      "per_primer": {
        "<primer_name>": {
          "found": int, "kept": int,
          "filtered_too_small": int, "filtered_too_large": int,
          "contigs": ["NC_xxx", ...]
        }
      },
      "per_contig": {
        "<sequence_id>": {"found": int, "kept": int}
      },
      "amplicons": [<amplicon_record>, ...]   # kept only
    }

Each <amplicon_record> contains:
    primer, amplimer_number, sequence_id, sequence_desc, length,
    forward: {primer_seq, strand, pos, mismatches},
    reverse: {primer_seq, strand, pos, mismatches},
    filter_status: "kept" | "too_small" | "too_large"
"""

import argparse
import json
import re
import sys

PRIMER_HIT_RE = re.compile(
    r'^\s*(?P<seq>\S+)\s+hits\s+(?P<strand>forward|reverse)\s+strand\s+at\s+'
    r'\[?(?P<pos>\d+)\]?\s+with\s+(?P<mm>\d+)\s+mismatches\s*$'
)
AMPLIMER_LENGTH_RE = re.compile(r'Amplimer length:\s+(\d+)\s+bp')
AMPLIMER_HEADER_RE = re.compile(r'^Amplimer\s+(\d+)\s*$')


def parse_primersearch_output(input_file):
    amplicons = []
    current_primer = None
    current = None

    def flush():
        nonlocal current
        if current and current.get('length') is not None:
            amplicons.append(current)
        current = None

    with open(input_file, 'r') as f:
        prev_was_sequence = False
        for raw in f:
            line = raw.rstrip('\n')
            stripped = line.strip()

            if line.startswith('Primer name '):
                flush()
                name = line[len('Primer name '):].strip()
                current_primer = name if name else None
                continue

            m = AMPLIMER_HEADER_RE.match(stripped)
            if m:
                flush()
                current = {
                    'primer': current_primer,
                    'amplimer_number': int(m.group(1)),
                    'sequence_id': None,
                    'sequence_desc': None,
                    'length': None,
                    'forward': None,
                    'reverse': None,
                }
                prev_was_sequence = False
                continue

            if current is None:
                continue

            if stripped.startswith('Sequence:'):
                current['sequence_id'] = stripped[len('Sequence:'):].strip() or None
                prev_was_sequence = True
                continue

            if prev_was_sequence and stripped and not stripped.startswith('Amplimer') \
                    and 'hits forward strand' not in stripped \
                    and 'hits reverse strand' not in stripped \
                    and 'Amplimer length:' not in stripped:
                current['sequence_desc'] = stripped
                prev_was_sequence = False
                continue
            prev_was_sequence = False

            hit = PRIMER_HIT_RE.match(line)
            if hit:
                record = {
                    'primer_seq': hit.group('seq'),
                    'strand': hit.group('strand'),
                    'pos': int(hit.group('pos')),
                    'mismatches': int(hit.group('mm')),
                }
                # primersearch always reports the forward-strand hit first, then
                # the reverse-strand hit. Use strand as the slot key.
                if hit.group('strand') == 'forward':
                    current['forward'] = record
                else:
                    current['reverse'] = record
                continue

            lm = AMPLIMER_LENGTH_RE.search(line)
            if lm:
                current['length'] = int(lm.group(1))
                continue

        flush()

    return amplicons


def classify(amp, min_size, max_size):
    length = amp['length']
    if length < min_size:
        return 'too_small'
    if length > max_size:
        return 'too_large'
    return 'kept'


def build_stats(amplicons, min_size, max_size):
    per_primer = {}
    per_contig = {}
    totals = {
        'total_amplicons': len(amplicons),
        'retained': 0,
        'filtered_too_small': 0,
        'filtered_too_large': 0,
    }
    kept_records = []

    for amp in amplicons:
        primer = amp['primer'] or '_unnamed_'
        contig = amp['sequence_id'] or '_unknown_'
        status = classify(amp, min_size, max_size)
        amp['filter_status'] = status

        pp = per_primer.setdefault(primer, {
            'found': 0, 'kept': 0,
            'filtered_too_small': 0, 'filtered_too_large': 0,
            'contigs': set(),
        })
        pp['found'] += 1
        pp['contigs'].add(contig)

        pc = per_contig.setdefault(contig, {'found': 0, 'kept': 0})
        pc['found'] += 1

        if status == 'kept':
            totals['retained'] += 1
            pp['kept'] += 1
            pc['kept'] += 1
            kept_records.append(amp)
        elif status == 'too_small':
            totals['filtered_too_small'] += 1
            pp['filtered_too_small'] += 1
        else:
            totals['filtered_too_large'] += 1
            pp['filtered_too_large'] += 1

    for v in per_primer.values():
        v['contigs'] = sorted(v['contigs'])

    return {
        'min_size': min_size,
        'max_size': max_size,
        'totals': totals,
        'per_primer': per_primer,
        'per_contig': per_contig,
        'amplicons': kept_records,
    }


def write_filtered_output(kept_amplicons, output_file):
    """Re-emit kept amplicons in primersearch format for backward compatibility."""
    current_primer = None
    with open(output_file, 'w') as f:
        # group by primer
        by_primer = {}
        for amp in kept_amplicons:
            by_primer.setdefault(amp['primer'], []).append(amp)

        for primer, amps in by_primer.items():
            f.write(f"\nPrimer name {primer}\n")
            for amp in amps:
                f.write(f"Amplimer {amp['amplimer_number']}\n")
                if amp['sequence_id']:
                    f.write(f"\tSequence: {amp['sequence_id']}  \n")
                if amp['sequence_desc']:
                    f.write(f"\t{amp['sequence_desc']}\n")
                if amp['forward']:
                    fw = amp['forward']
                    f.write(f"\t{fw['primer_seq']} hits forward strand at {fw['pos']} with {fw['mismatches']} mismatches\n")
                if amp['reverse']:
                    rv = amp['reverse']
                    f.write(f"\t{rv['primer_seq']} hits reverse strand at [{rv['pos']}] with {rv['mismatches']} mismatches\n")
                f.write(f"\tAmplimer length: {amp['length']} bp\n")


def main():
    parser = argparse.ArgumentParser(description='Filter amplicons from primersearch output')
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--stats', required=True)
    parser.add_argument('--min-size', type=int, default=50)
    parser.add_argument('--max-size', type=int, default=10000)
    args = parser.parse_args()

    amplicons = parse_primersearch_output(args.input)
    stats = build_stats(amplicons, args.min_size, args.max_size)

    kept = [a for a in stats['amplicons']]
    write_filtered_output(kept, args.output)

    with open(args.stats, 'w') as f:
        json.dump(stats, f, indent=2)

    t = stats['totals']
    print(f"Filtered {t['total_amplicons']} amplicons:", file=sys.stderr)
    print(f"  - Retained: {t['retained']}", file=sys.stderr)
    print(f"  - Too small (<{args.min_size} bp): {t['filtered_too_small']}", file=sys.stderr)
    print(f"  - Too large (>{args.max_size} bp): {t['filtered_too_large']}", file=sys.stderr)
    print(f"  - Primers: {len(stats['per_primer'])}; contigs: {len(stats['per_contig'])}", file=sys.stderr)


if __name__ == '__main__':
    main()
