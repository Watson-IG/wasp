#!/usr/bin/env python3
import os
import sys
import pandas as pd
import glob

def get_locus_from_gene(gene):
    for loc in ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]:
        if gene.startswith(loc):
            return loc
    return None

def main():
    if len(sys.argv) < 3:
        print("Usage: python check_missing_genes.py <alleles_dir> <gene.bed>")
        sys.exit(1)
        
    alleles_dir = sys.argv[1]
    gene_bed = sys.argv[2]
    
    if not os.path.exists(gene_bed):
        # Silently exit if gene.bed is not found
        sys.exit(0)
        
    # Read gene.bed and group expected genes by locus
    expected_genes_by_locus = {}
    with open(gene_bed, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                gene = parts[3]
                locus = get_locus_from_gene(gene)
                if locus:
                    if locus not in expected_genes_by_locus:
                        expected_genes_by_locus[locus] = set()
                    expected_genes_by_locus[locus].add(gene)
                    
    # Group available allele CSV files by locus
    all_files = glob.glob(os.path.join(alleles_dir, "*.csv"))
    locus_files = {}
    for f in all_files:
        basename = os.path.basename(f)
        for loc in ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]:
            if f"_{loc}_" in basename:
                if loc not in locus_files:
                    locus_files[loc] = []
                locus_files[loc].append(f)
                break
                
    # Check each locus for missing genes
    for locus, genes in expected_genes_by_locus.items():
        if locus not in locus_files:
            continue
            
        # Determine the final table: combined > reference > digger
        files = locus_files[locus]
        final_file = None
        for ftype in ["combined", "reference", "digger"]:
            for f in files:
                if f"_{ftype}_" in os.path.basename(f):
                    final_file = f
                    break
            if final_file:
                break
                
        if not final_file:
            continue
            
        try:
            df = pd.read_csv(final_file)
            if 'gene' in df.columns:
                found_genes = set(df['gene'].dropna().astype(str))
                missing_genes = genes - found_genes
                for mg in sorted(missing_genes):
                    print(f"WARNING: Missing genes for: {mg}")
        except Exception as e:
            pass

if __name__ == "__main__":
    main()
