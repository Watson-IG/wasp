import os
import argparse
import pandas as pd
import glob
import gzip
import concurrent.futures

def process_file(directory, pattern, fasta_output, name_col, seq_col, filter_fn):
    pattern_path = os.path.join(directory, pattern)
    matching_files = glob.glob(pattern_path)
    
    if not matching_files:
        print(f"Warning: No files found for pattern {pattern}. Skipping.")
        return
    
    input_path = matching_files[0]
    if len(matching_files) > 1:
        print(f"Multiple files found for pattern {pattern}. Using {os.path.basename(input_path)}")
    
    fasta_output_path = os.path.join(directory, fasta_output)
    
    # Read CSV file
    try:
        df = pd.read_csv(input_path)
    except Exception as e:
        print(f"Error reading {input_path}: {e}")
        return
    
    # Check required columns exist
    if name_col not in df.columns or seq_col not in df.columns:
        print(f"Warning: Required columns ('{name_col}', '{seq_col}') not found in {input_path}. Skipping.")
        return
    
    # Apply filtering
    filtered = filter_fn(df)
    
    # Remove duplicates
    unique_entries = filtered.drop_duplicates(subset=[name_col, seq_col])
    
    # Write FASTA file as gzip
    try:
        with gzip.open(fasta_output_path, "wt") as fasta:
            for _, row in unique_entries.iterrows():
                fasta.write(f">{row[name_col]}\n{row[seq_col]}\n")
    except Exception as e:
        print(f"Error writing {fasta_output_path}: {e}")
        return
    
    print(f"Processed {os.path.basename(input_path)} -> {fasta_output}")


def ref_filter(df):
    """Filter for reference-guided and combined allele tables."""
    is_c_gene = df["vdjbase_allele"].str[3] == "C"
    c_filter = df["Average_Coverage"] >= 30
    non_c_filter = (
        (df["Average_Coverage"] >= 30) &
        (df["Fully_Spanning_Reads_100%_Match"] >= 10)
    )
    combined_filter = (is_c_gene & c_filter) | (~is_c_gene & non_c_filter)
    return df[combined_filter]


def digger_filter(df):
    """Filter for digger allele tables (no coverage filtering for now)."""
    return df.dropna(subset=["gene type", "seq"])


# Mode-specific configuration: (file pattern suffix, output suffix, name_col, seq_col, filter_fn)
MODE_CONFIG = {
    "ref": {
        "name_col": "vdjbase_allele",
        "seq_col": "gene_sequence",
        "filter_fn": ref_filter,
        "loci": [
            ("*_IGHC_reference_annotated-alleles.csv", "IGHC_alleles.fasta.gz"),
            ("*_TRA_reference_annotated-alleles.csv",  "TRA_alleles.fasta.gz"),
            ("*_IGH_reference_annotated-alleles.csv",  "IGH_alleles.fasta.gz"),
            ("*_TRB_reference_annotated-alleles.csv",  "TRB_alleles.fasta.gz"),
            ("*_IGK_reference_annotated-alleles.csv",  "IGK_alleles.fasta.gz"),
            ("*_TRD_reference_annotated-alleles.csv",  "TRD_alleles.fasta.gz"),
            ("*_IGL_reference_annotated-alleles.csv",  "IGL_alleles.fasta.gz"),
            ("*_TRG_reference_annotated-alleles.csv",  "TRG_alleles.fasta.gz"),
        ],
    },
    "digger": {
        "name_col": "gene type",
        "seq_col": "seq",
        "filter_fn": digger_filter,
        "loci": [
            ("*_IGH_digger_annotated-alleles.csv", "IGH_alleles.fasta.gz"),
            ("*_IGK_digger_annotated-alleles.csv", "IGK_alleles.fasta.gz"),
            ("*_IGL_digger_annotated-alleles.csv", "IGL_alleles.fasta.gz"),
            ("*_TRA_digger_annotated-alleles.csv", "TRA_alleles.fasta.gz"),
            ("*_TRB_digger_annotated-alleles.csv", "TRB_alleles.fasta.gz"),
            ("*_TRD_digger_annotated-alleles.csv", "TRD_alleles.fasta.gz"),
            ("*_TRG_digger_annotated-alleles.csv", "TRG_alleles.fasta.gz"),
        ],
    },
    "combined": {
        "name_col": "vdjbase_allele",
        "seq_col": "gene_sequence",
        "filter_fn": ref_filter,
        "loci": [
            ("*_IGH_combined_annotated-alleles.csv", "IGH_alleles.fasta.gz"),
            ("*_IGK_combined_annotated-alleles.csv", "IGK_alleles.fasta.gz"),
            ("*_IGL_combined_annotated-alleles.csv", "IGL_alleles.fasta.gz"),
            ("*_TRA_combined_annotated-alleles.csv", "TRA_alleles.fasta.gz"),
            ("*_TRB_combined_annotated-alleles.csv", "TRB_alleles.fasta.gz"),
            ("*_TRD_combined_annotated-alleles.csv", "TRD_alleles.fasta.gz"),
            ("*_TRG_combined_annotated-alleles.csv", "TRG_alleles.fasta.gz"),
        ],
    },
}


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Convert CSV allele files to FASTA format using multithreading.'
    )
    parser.add_argument(
        '--directory', 
        type=str, 
        default='.', 
        help='Directory containing input CSV files (default: current directory)'
    )
    parser.add_argument(
        '--mode',
        type=str,
        choices=['ref', 'digger', 'combined', 'denovo'],
        default='ref',
        help='Pipeline mode: ref, digger/denovo, or combined (default: ref)'
    )
    args = parser.parse_args()

    # Map 'denovo' to 'digger' since they use the same files
    mode = 'digger' if args.mode == 'denovo' else args.mode

    config = MODE_CONFIG[mode]
    name_col = config["name_col"]
    seq_col = config["seq_col"]
    filter_fn = config["filter_fn"]
    files_info = config["loci"]
    
    # Use ThreadPoolExecutor to process files concurrently
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for pattern, fasta_output in files_info:
            futures.append(
                executor.submit(process_file, args.directory, pattern, fasta_output, name_col, seq_col, filter_fn)
            )
        # Wait for all threads to complete
        concurrent.futures.wait(futures)

if __name__ == "__main__":
    main()
