#!/usr/bin/env python3
import argparse
import os
import subprocess
import tempfile
from collections import defaultdict
import pysam

def main():
    parser = argparse.ArgumentParser(description="Cut contigs using reference-independent BLAST against VDJ alleles.")
    parser.add_argument("--fasta", required=True, help="Input assembly FASTA")
    parser.add_argument("--allele_ref_dir", required=True, help="Directory containing VDJ reference FASTAs")
    parser.add_argument("--buffer", type=int, default=20000, help="Buffer distance to add to left and right flanks")
    parser.add_argument("--out", required=True, help="Output cut assembly FASTA")
    args = parser.parse_args()

    # Step 1: Create a combined reference FASTA
    temp_dir = tempfile.mkdtemp()
    combined_ref = os.path.join(temp_dir, "combined_vdj_refs.fasta")
    
    loci_fastas = []
    for f in os.listdir(args.allele_ref_dir):
        if f.endswith(".fasta") or f.endswith(".fa"):
            loci_fastas.append(os.path.join(args.allele_ref_dir, f))
            
    if not loci_fastas:
        print(f"Error: No FASTA files found in {args.allele_ref_dir}")
        return

    with open(combined_ref, "w") as out_f:
        for fa in loci_fastas:
            with open(fa, "r") as in_f:
                out_f.write(in_f.read())

    # Step 2: Make BLAST DB
    db_path = os.path.join(temp_dir, "vdj_db")
    subprocess.run(
        ["makeblastdb", "-in", combined_ref, "-dbtype", "nucl", "-out", db_path],
        check=True, stdout=subprocess.DEVNULL
    )

    # Step 3: BLAST contigs against the DB
    blast_out = os.path.join(temp_dir, "blast_results.tsv")
    subprocess.run(
        [
            "blastn",
            "-query", args.fasta,
            "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore qstart qend slen",
            "-evalue", "1e-5",
            "-out", blast_out,
        ],
        check=True
    )

    # Step 4: Parse hits to find leftmost and rightmost bounds per contig+locus
    # contig -> locus -> [min_coord, max_coord]
    bounds = defaultdict(lambda: defaultdict(lambda: [float('inf'), float('-inf')]))
    
    known_loci = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]

    with open(blast_out) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            qseqid = fields[0]
            sseqid = fields[1]
            qstart = int(fields[6])
            qend = int(fields[7])
            slen = int(fields[8])
            
            # Enforce 95% of subject sequence length
            pident = float(fields[2])
            length = int(fields[3])
            if length < 0.95 * slen or pident < 80:
                continue

            # Determine locus
            assigned_locus = None
            for loc in known_loci:
                if loc in sseqid:
                    assigned_locus = loc
                    break
            
            if assigned_locus:
                # Group TRA and TRD together to prevent overlapping duplicate cuts
                if assigned_locus in ["TRA", "TRD"]:
                    assigned_locus = "TRA_TRD"
                    
                actual_start = min(qstart, qend)
                actual_end = max(qstart, qend)
                
                if actual_start < bounds[qseqid][assigned_locus][0]:
                    bounds[qseqid][assigned_locus][0] = actual_start
                if actual_end > bounds[qseqid][assigned_locus][1]:
                    bounds[qseqid][assigned_locus][1] = actual_end

    # Step 5: Cut contigs and write out
    fa = pysam.FastaFile(args.fasta)
    
    with open(args.out, "w") as out_f:
        for qseqid in bounds:
            if qseqid not in fa.references:
                continue
                
            seq = fa.fetch(qseqid)
            seq_len = len(seq)
            
            for locus, (leftmost, rightmost) in bounds[qseqid].items():
                cut_start = max(0, leftmost - args.buffer - 1) # 0-indexed
                cut_end = min(seq_len, rightmost + args.buffer)
                
                # Suffix the contig name with locus if there are multiple loci on this contig,
                # or always suffix it for clarity. Let's always suffix it.
                new_qseqid = f"{qseqid}_{locus}"
                
                cut_seq = seq[cut_start:cut_end]
                out_f.write(f">{new_qseqid}\n{cut_seq}\n")
                
    print(f"Cut {len(bounds)} contigs into {sum(len(l) for l in bounds.values())} sub-contigs.")

if __name__ == "__main__":
    main()
