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
    'v_intron': 'V-INTRON',
    'v_spacer': 'V-SPACER',
    'd_3_heptamer': 'D-3_HEPTAMER',
    'd_3_nonamer': 'D-3_NONAMER',
    'd_3_spacer': 'D-3_SPACER',
    'd_5_heptamer': 'D-5_HEPTAMER',
    'd_5_nonamer': 'D-5_NONAMER',
    'd_5_spacer': 'D-5_SPACER',
    'j_heptamer': 'J-HEPTAMER',
    'j_nonamer': 'J-NONAMER',
    'j_spacer': 'J-SPACER',
    'l_part1': 'L-PART1',
    'l_part2': 'L-PART2',
    'start': 'REGION_start',
    'end': 'REGION_end',
    'ref_match': 'closest_refdb_allele',
    'gene_seq': 'gene_sequence',
}


def merge_tables(ref_csv: str, digger_csv: str, output_csv: str) -> None:
    """Read the two source CSVs, tag with 'source', rename digger columns, and merge."""

    df_ref = pd.read_csv(ref_csv)
    df_ref['source'] = 'reference-guided'
    
    if 'closest_imgt_allele' in df_ref.columns:
        df_ref = df_ref.rename(columns={'closest_imgt_allele': 'closest_refdb_allele'})

    if digger_csv.endswith('.gz'):
        df_digger = pd.read_csv(digger_csv, compression='gzip')
    else:
        df_digger = pd.read_csv(digger_csv)
    df_digger['source'] = 'digger'
    
    if 'ref_match' in df_digger.columns:
        df_digger['gene'] = df_digger['ref_match'].str.split('*').str[0]

    if 'gene' in df_digger.columns:
        if 'seq' in df_digger.columns:
            is_v = df_digger['gene'].str.match(r'^[A-Z]{3}V', na=False)
            is_d = df_digger['gene'].str.match(r'^[A-Z]{3}D', na=False)
            is_j = df_digger['gene'].str.match(r'^[A-Z]{3}J', na=False)
            df_digger.loc[is_v, 'V-REGION'] = df_digger.loc[is_v, 'seq']
            df_digger.loc[is_d, 'D-REGION'] = df_digger.loc[is_d, 'seq']
            df_digger.loc[is_j, 'J-REGION'] = df_digger.loc[is_j, 'seq']
            df_digger = df_digger.drop(columns=['seq'])
            
        if 'seq_gapped' in df_digger.columns:
            is_v = df_digger['gene'].str.match(r'^[A-Z]{3}V', na=False)
            df_digger.loc[is_v, 'V-REGION-GAPPED'] = df_digger.loc[is_v, 'seq_gapped']
            df_digger = df_digger.drop(columns=['seq_gapped'])

    df_digger = df_digger.rename(columns=DIGGER_TO_REF_RENAME)
    
    for col in ['project', 'subject', 'sample_name']:
        if col in df_ref.columns:
            val = df_ref[col].dropna().iloc[0] if not df_ref[col].dropna().empty else None
            df_digger[col] = val

    if 'contig' in df_digger.columns:
        df_digger['haplotype'] = df_digger['contig'].astype(str).str.extract(r'_hap(\d+)', expand=False)
    
    if 'V-REGION' in df_ref.columns and 'V-REGION' in df_digger.columns and 'contig' in df_ref.columns and 'contig' in df_digger.columns:
        ref_valid = df_ref.dropna(subset=['V-REGION', 'contig'])
        ref_pairs = set(zip(ref_valid['V-REGION'], ref_valid['contig']))
        
        mask = df_digger.apply(lambda row: (row['V-REGION'], row['contig']) in ref_pairs, axis=1)
        df_digger = df_digger[~mask]
    elif 'V-REGION' in df_ref.columns and 'V-REGION' in df_digger.columns:
        ref_v_regions = set(df_ref['V-REGION'].dropna())
        df_digger = df_digger[~df_digger['V-REGION'].isin(ref_v_regions)]

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
