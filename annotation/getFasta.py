import os
import argparse
import pandas as pd
import glob
import gzip
import concurrent.futures

def process_file(directory, pattern, fasta_output):
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
    
    # Apply filtering
    is_c_gene = df["vdjbase_allele"].str[3] == "C"
    c_filter = df["Average_Coverage"] >= 30
    non_c_filter = (
        (df["Average_Coverage"] >= 30) &
        (df["Fully_Spanning_Reads_100%_Match"] >= 10)
    )
    combined_filter = (is_c_gene & c_filter) | (~is_c_gene & non_c_filter)
    filtered = df[combined_filter]
    
    # Remove duplicates
    unique_entries = filtered.drop_duplicates(subset=["vdjbase_allele", "gene_sequence"])
    
    # Write FASTA file as gzip
    try:
        with gzip.open(fasta_output_path, "wt") as fasta:
            for _, row in unique_entries.iterrows():
                fasta.write(f">{row['vdjbase_allele']}\n{row['gene_sequence']}\n")
    except Exception as e:
        print(f"Error writing {fasta_output_path}: {e}")
        return
    
    print(f"Processed {os.path.basename(input_path)} -> {fasta_output}")

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
    args = parser.parse_args()
    
    # Define input file patterns and corresponding gzip FASTA output files
    files_info = [
        ("*_IGHC_annotated-alles-with-read-support.csv", "IGHC_alleles.fasta.gz"),
        ("*_TRA_annotated-alles-with-read-support.csv", "TRA_alleles.fasta.gz"),
        ("*_IGH_annotated-alles-with-read-support.csv", "IGH_alleles.fasta.gz"),
        ("*_TRB_annotated-alles-with-read-support.csv", "TRB_alleles.fasta.gz"),
        ("*_IGK_annotated-alles-with-read-support.csv", "IGK_alleles.fasta.gz"),
        ("*_TRD_annotated-alles-with-read-support.csv", "TRD_alleles.fasta.gz"),
        ("*_IGL_annotated-alles-with-read-support.csv", "IGL_alleles.fasta.gz"),
        ("*_TRG_annotated-alles-with-read-support.csv", "TRG_alleles.fasta.gz")
    ]
    
    # Use ThreadPoolExecutor to process files concurrently
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for pattern, fasta_output in files_info:
            futures.append(
                executor.submit(process_file, args.directory, pattern, fasta_output)
            )
        # Wait for all threads to complete
        concurrent.futures.wait(futures)

if __name__ == "__main__":
    main()
