#!/usr/bin/env python3
"""
merge_allele_tables.py

Merges reference-guided and digger allele annotation tables per locus for
combined mode. Produces a merged CSV with a 'source' column indicating the
origin of each row.

Usage:
    python merge_allele_tables.py <ref_guided_csv> <digger_csv> <output_csv>
"""

import argparse
import os
import sys
import pandas as pd


# Map digger column names to reference-guided column names so that
# the two tables can be concatenated cleanly.
DIGGER_TO_REF_RENAME = {
    'v_heptamer': 'V-HEPTAMER',
    'v_nonamer': 'V-NONAMER',
    'd_3_heptamer': 'D-3_HEPTAMER',
    'd_3_nonamer': 'D-3_NONAMER',
    'd_5_heptamer': 'D-5_HEPTAMER',
    'd_5_nonamer': 'D-5_NONAMER',
    'j_heptamer': 'J-HEPTAMER',
    'j_nonamer': 'J-NONAMER',
    'l_part1': 'L-PART1',
    'l_part2': 'L-PART2',
    'seq': 'V-REGION',
    'seq_gapped': 'V-REGION-GAPPED',
    'start': 'REGION_start',
    'end': 'REGION_end',
}


def merge_tables(ref_csv: str, digger_csv: str, output_csv: str) -> None:
    """Read the two source CSVs, tag with 'source', rename digger columns, and merge."""

    df_ref = pd.read_csv(ref_csv)
    df_ref['source'] = 'reference-guided'

    if digger_csv.endswith('.gz'):
        df_digger = pd.read_csv(digger_csv, compression='gzip')
    else:
        df_digger = pd.read_csv(digger_csv)
    df_digger['source'] = 'digger'

    df_digger = df_digger.rename(columns=DIGGER_TO_REF_RENAME)

    merged = pd.concat([df_ref, df_digger], axis=0, ignore_index=True)
    merged.to_csv(output_csv, index=False)
    print(f"Merged {len(df_ref)} reference-guided + {len(df_digger)} digger rows -> {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description='Merge reference-guided and digger allele tables for combined mode.'
    )
    parser.add_argument('ref_csv', help='Path to the reference-guided allele CSV')
    parser.add_argument('digger_csv', help='Path to the digger allele CSV')
    parser.add_argument('output_csv', help='Path for the merged output CSV')
    args = parser.parse_args()

    for path, label in [(args.ref_csv, 'Reference-guided CSV'), (args.digger_csv, 'Digger CSV')]:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {label} not found: {path}")

    merge_tables(args.ref_csv, args.digger_csv, args.output_csv)


if __name__ == '__main__':
    main()
