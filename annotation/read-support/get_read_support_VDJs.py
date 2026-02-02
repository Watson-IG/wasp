#!/usr/bin/env python3
"""
This script reproduces the bioinformatics pipeline logic in a single Python file.
It integrates previously external scripts and includes robust sanity checks.

Dependencies:
    - pandas
    - pysam
    - bedtools (system binary)
    - samtools (system binary)
    - minimap2 (system binary)
"""

import sys
import os
import subprocess
import csv
import shutil
import pandas as pd
import pysam

# ==========================================
# Helper Functions & Sanity Checks
# ==========================================

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in seq[::-1]])

def validate_vdj_stats(sample_name, gene_name, stats_tuple):
    """
    Sanity Check: Validates that the summary statistics (columns) match the 
    raw colon-separated position data strings.
    
    Args:
        sample_name (str): The name of the sample.
        gene_name (str): The gene being analyzed.
        stats_tuple (tuple): Output from parse_mpileup_and_calculate containing:
            0: Total_Positions (int)
            1: Avg_Coverage (float) - Unused here
            2: Mismatched_Positions (int) - Count of positions with >20% mismatch
            3: Matched_Positions (int)    - Count of positions with >80% match
            4: Position_Mismatches (str)  - Colon-separated string (e.g., "0:1:0")
            5: Position_Matches (str)     - Colon-separated string (e.g., "10:9:10")
            6: Percent_Accuracy (float)   - Unused here
            7: Positions_With_10x (int)   - Count of positions with coverage >= 10
    """
    
    # --- Step 1: Unpack Reported Stats from the Tuple ---
    expected_total_positions = stats_tuple[0]
    reported_mismatch_count  = stats_tuple[2]
    reported_match_count     = stats_tuple[3]
    reported_10x_count       = stats_tuple[7]
    
    raw_mismatch_string = stats_tuple[4]
    raw_match_string    = stats_tuple[5]
    
    # Convert colon-separated strings into lists of numbers
    # Example: "10:12:0" -> ['10', '12', '0']
    mismatch_list = raw_mismatch_string.split(':') if raw_mismatch_string else []
    match_list    = raw_match_string.split(':')    if raw_match_string    else []

    # --- Step 2: verify Data Length Consistency ---
    # The number of data points in the string must match the total region length
    if len(mismatch_list) != expected_total_positions or len(match_list) != expected_total_positions:
        print(f"[SANITY FAIL] {sample_name} {gene_name}: Data gap detected!")
        print(f"  Expected {expected_total_positions} positions based on coordinates.")
        print(f"  Found {len(match_list)} match entries and {len(mismatch_list)} mismatch entries.")
        return False

    # --- Step 3: Recalculate Stats from Raw Data ---
    calculated_match_count = 0
    calculated_mismatch_count = 0
    calculated_10x_count = 0

    # Iterate through every position in the region
    for match_str_val, mismatch_str_val in zip(match_list, mismatch_list):
        try:
            # Parse read counts (int(float()) handles cases like "10.0" safely)
            match_read_count = int(float(match_str_val))
            mismatch_read_count = int(float(mismatch_str_val))
            
            total_coverage_at_pos = match_read_count + mismatch_read_count
            
            if total_coverage_at_pos > 0:
                # Check definition: Matches must be > 80% of reads
                match_ratio = match_read_count / total_coverage_at_pos
                if match_ratio > 0.8:
                    calculated_match_count += 1
                
                # Check definition: Mismatches must be > 20% of reads
                mismatch_ratio = mismatch_read_count / total_coverage_at_pos
                if mismatch_ratio > 0.2:
                    calculated_mismatch_count += 1
            
            # Check definition: Coverage must be >= 10 reads
            if total_coverage_at_pos >= 10:
                calculated_10x_count += 1
                
        except ValueError:
            # Skip positions with malformed data
            continue

    # --- Step 4: Compare Reported vs Calculated ---
    errors = []
    
    if reported_match_count != calculated_match_count:
        errors.append(f"Matched_Positions mismatch: Header says {reported_match_count}, Raw Data says {calculated_match_count}")
        
    if reported_mismatch_count != calculated_mismatch_count:
        errors.append(f"Mismatched_Positions mismatch: Header says {reported_mismatch_count}, Raw Data says {calculated_mismatch_count}")
        
    if reported_10x_count != calculated_10x_count:
        errors.append(f"10x_Coverage mismatch: Header says {reported_10x_count}, Raw Data says {calculated_10x_count}")

    # --- Step 5: Report Errors ---
    if errors:
        print(f"[SANITY FAIL] {sample_name} {gene_name} Logic Error:")
        for error in errors:
            print(f"  - {error}")
        return False

    return True

def append_pos_import_genes_internal(input_csv, original_output_csv, modified_output_csv):
    """Integrates logic from append_pos_import_genes.py"""
    shutil.copyfile(original_output_csv, modified_output_csv)
    try:
        output_df = pd.read_csv(modified_output_csv)
        if 'genotyper_gene' in output_df.columns:
            output_df.rename(columns={'genotyper_gene': 'gene'}, inplace=True)

        input_df = pd.read_csv(input_csv)
        for _, input_row in input_df.iterrows():
            gene_name = input_row['gene']
            if any(s in gene_name for s in ['IGKV', 'IGLV', 'IGHV', 'TRBV', 'TRDV', 'TRGV', 'TRAV']):
                start_col, end_col = 'V-REGION_start', 'V-REGION_end'
            elif any(s in gene_name for s in ['IGKJ', 'IGLJ', 'IGHJ', 'TRBJ', 'TRDJ', 'TRGJ', 'TRAJ']):
                start_col, end_col = 'J-REGION_start', 'J-REGION_end'
            elif any(s in gene_name for s in ['IGHD', 'TRBD', 'TRDD']):
                start_col, end_col = 'D-REGION_start', 'D-REGION_end'
            elif any(s in gene_name for s in ['IGKC', 'IGLC', 'TRBC', 'TRDC', 'TRGC', 'TRAC']):
                start_col, end_col = 'allele_sequence_start', 'allele_sequence_end'
            else:
                continue

            if pd.notna(input_row.get(start_col)) and pd.notna(input_row.get(end_col)):
                mask = (output_df['contig'] == input_row['contig']) & (output_df['gene'] == input_row['gene'])
                if mask.any():
                    output_df.loc[mask, 'REGION_start'] = input_row[start_col]
                    output_df.loc[mask, 'REGION_end'] = input_row[end_col]

        if 'notes' in output_df.columns:
            output_df['notes'] = output_df['notes'].fillna('').astype(str).str.replace(',', ';', regex=False)
        output_df.to_csv(modified_output_csv, index=False)
    except Exception as e:
        print(f"[WARNING] Error in append_pos_import_genes_internal: {e}")

def ighc_append_pos_internal(input_csv, original_output_csv, modified_output_csv):
    """Integrates logic from ighc_append_pos.py"""
    shutil.copyfile(original_output_csv, modified_output_csv)
    try:
        output_df = pd.read_csv(modified_output_csv)
        if 'genotyper_gene' in output_df.columns:
            output_df.rename(columns={'genotyper_gene': 'gene'}, inplace=True)
        input_df = pd.read_csv(input_csv)
        
        exon_columns = [f'C-EXON_{i}_{s}' for i in range(1, 10) for s in ['start', 'end']]

        for _, input_row in input_df.iterrows():
            mask = (output_df['contig'] == input_row['contig']) & (output_df['gene'] == input_row['gene'])
            if not mask.any(): continue

            min_start = float('inf')
            max_end = 0

            for col in exon_columns:
                if col in input_row and pd.notna(input_row[col]):
                    val = input_row[col]
                    output_df.loc[mask, col] = val
                    if col.endswith("_start") and val > 0:
                        min_start = min(min_start, val)
                    elif col.endswith("_end"):
                        max_end = max(max_end, val)

            if min_start < float('inf'):
                output_df.loc[mask, 'allele_sequence_start'] = min_start
            if max_end > 0:
                output_df.loc[mask, 'allele_sequence_end'] = max_end

        output_df.to_csv(modified_output_csv, index=False)
    except Exception as e:
        print(f"[WARNING] Error in ighc_append_pos_internal: {e}")

def get_sequence_from_row(row, gene_key):
    """Extracts sequence string from CSV row."""
    reverse_comp = False
    if 'sense' in row and row['sense'] == '-':
        reverse_comp = True

    sequence_column = None
    if any(s in gene_key for s in ['IGKV', 'IGLV', 'IGHV', 'TRAV', 'TRBV', 'TRDV', 'TRGV']):
        sequence_column = 'V-REGION'
    elif any(s in gene_key for s in ['IGKJ', 'IGLJ', 'IGHJ', 'TRAJ', 'TRBJ', 'TRDJ', 'TRGJ']):
        sequence_column = 'J-REGION'
    elif any(s in gene_key for s in ['IGHD', 'TRBD', 'TRDD']):
        sequence_column = 'D-REGION'
    elif any(s in gene_key for s in ['IGKC', 'IGLC', 'TRAC', 'TRBC', 'TRDC', 'TRGC']):
        sequence_column = 'C-REGION'
    
    if not sequence_column or sequence_column not in row:
        return ""

    sequence = row[sequence_column]
    if reverse_comp and sequence:
        sequence = reverse_complement(sequence)
    return sequence

def count_matching_reads_internal(bamfile, chrom, start, end, sequence):
    """Counts full span and perfect matches using pysam."""
    full_span_count = 0
    perfect_match_count = 0
    try:
        samfile = pysam.AlignmentFile(bamfile, "rb")
        for read in samfile.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_start <= start and read.reference_end >= end:
                full_span_count += 1
                read_seq = read.query_sequence
                if sequence and (sequence in read_seq or reverse_complement(sequence) in read_seq):
                    perfect_match_count += 1
        samfile.close()
    except Exception as e:
        print(f"[WARNING] Pysam error in count_matching_reads: {e}")
        return 0, 0
    return full_span_count, perfect_match_count

def get_ighc_sequences_from_row(row):
    """Extracts IGHC exon sequences."""
    sequences = []
    exon_start_end_list = []
    reverse_comp = 'sense' in row and row['sense'] == '-'
    
    for i in range(1, 10):
        seq_col = f'C-EXON_{i}'
        start_col = f'C-EXON_{i}_start'
        end_col = f'C-EXON_{i}_end'
        sequence = row.get(seq_col, '')
        start_val = row.get(start_col, '')
        end_val = row.get(end_col, '')

        if start_val and end_val:
            try:
                s_int = int(float(start_val))
                e_int = int(float(end_val))
                exon_start_end_list.append((s_int, e_int))
            except ValueError: pass
        
        if reverse_comp and sequence:
            sequence = reverse_complement(sequence)
        if sequence:
            sequences.append(sequence)
            
    return sequences, exon_start_end_list

def count_ighc_matches_internal(bamfile, contig, start, end, sequences, exon_start_end_list):
    """IGHC matching logic."""
    full_span_count = 0
    full_span_all_match_count = 0
    
    max_len = max(len(sequences), len(exon_start_end_list))
    if max_len == 0: return 0, 0, [], []
    
    perfect_matches = [0] * len(sequences)
    perfect_spans = [0] * len(exon_start_end_list)

    try:
        samfile = pysam.AlignmentFile(bamfile, "rb")
        for read in samfile.fetch(contig, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            
            if read.reference_start <= start and read.reference_end >= end:
                full_span_count += 1
                read_seq = read.query_sequence
                all_matches_flag = True
                for seq in sequences:
                    if not seq: continue
                    if (seq not in read_seq) and (reverse_complement(seq) not in read_seq):
                        all_matches_flag = False
                        break
                if all_matches_flag and sequences:
                    full_span_all_match_count += 1
            
            read_seq = read.query_sequence
            for idx, (ex_s, ex_e) in enumerate(exon_start_end_list):
                if read.reference_start <= ex_s and read.reference_end >= ex_e:
                    perfect_spans[idx] += 1
            for idx, seq in enumerate(sequences):
                if seq and (seq in read_seq or reverse_complement(seq) in read_seq):
                    perfect_matches[idx] += 1
        samfile.close()
    except Exception as e:
        print(f"[WARNING] Pysam error in ighc count: {e}")

    return full_span_count, full_span_all_match_count, perfect_matches, perfect_spans

# ==========================================
# Main Pipeline Logic
# ==========================================

def safe_run(cmd_args):
    cmd_args = [str(arg) for arg in cmd_args if str(arg).strip()]
    try:
        result = subprocess.run(cmd_args, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[WARNING] Command failed: {' '.join(cmd_args)}")
            print(f"stderr:\n{result.stderr}")
        return result
    except Exception as e:
        print(f"[WARNING] Exception: {' '.join(cmd_args)}")
        print(e)
    return None

def run_make_ref_masked(reffn, IG_loci, scratch):
    masked_ref = os.path.join(scratch, "ref_IG_masked.fasta")
    safe_run(["bedtools", "maskfasta", "-fi", reffn, "-bed", IG_loci, "-fo", masked_ref])
    safe_run(["samtools", "faidx", masked_ref])

def run_map_ccs_to_pers(fofn, scratch, mask_ref, minimap_option, threads):
    reads_fasta = os.path.join(scratch, "reads.fasta")
    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for cols in reader:
            if not cols or len(cols) <= 18: continue
            sample = cols[0]
            asm_bam = cols[1]
            print(f"Sample is {sample}")
            outd = os.path.join(scratch, "read_support", sample)
            os.makedirs(os.path.join(outd, "ccs_to_pers"), exist_ok=True)

            safe_run(["samtools", "faidx", reads_fasta])
            contigs_fa = os.path.join(outd, "ccs_to_pers", "contigs.fasta")
            try:
                with subprocess.Popen(["samtools", "view", "-F", "0x100", "-F", "0x800", asm_bam], stdout=subprocess.PIPE, text=True) as proc, open(contigs_fa, 'w') as contigs_out:
                    for line in proc.stdout:
                        fields = line.strip().split('\t')
                        if len(fields) >= 10: contigs_out.write(f">{fields[0]}\n{fields[9]}\n")
            except Exception as e: print(f"[WARNING] Contig extraction failed: {e}")
            
            safe_run(["samtools", "faidx", contigs_fa])
            pers_ref = os.path.join(outd, "ccs_to_pers", "pers_ref.fasta")
            try:
                with open(pers_ref, 'w') as p, open(mask_ref, 'r') as m, open(contigs_fa, 'r') as c:
                    p.write(m.read() + c.read())
            except Exception: pass
            safe_run(["samtools", "faidx", pers_ref])
            
            sam_out = os.path.join(outd, "ccs_to_pers", "output.sam")
            bam_out = os.path.join(outd, "ccs_to_pers", "output.bam")
            sorted_bam = os.path.join(outd, "ccs_to_pers", "output.sorted.bam")
            safe_run(["minimap2", "-ax", minimap_option, "--secondary=yes", "-t", str(threads), "-L", pers_ref, reads_fasta, "-o", sam_out])
            safe_run(["samtools", "view", "-Sbh", sam_out, "-o", bam_out])
            safe_run(["samtools", "sort", "-@", str(threads), bam_out, "-o", sorted_bam])
            safe_run(["samtools", "index", sorted_bam])
            try: os.remove(sam_out)
            except OSError: pass

def run_append_pos(fofn, scratch):
    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for cols in reader:
            if not cols or len(cols) < 18: continue
            sample = cols[0]
            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")
            
            # Map paths
            paths = {
                "IGK": {"src": cols[3], "out": os.path.join(base_outd, "IGK", os.path.basename(cols[3]))},
                "IGL": {"src": cols[5], "out": os.path.join(base_outd, "IGL", os.path.basename(cols[5]))},
                "IGH": {"src": cols[7], "out": os.path.join(base_outd, "IGH", os.path.basename(cols[7]))},
                "IGHC": {"src": cols[9], "out": os.path.join(base_outd, "IGHC", os.path.basename(cols[9]))},
                "TRB": {"src": cols[11], "out": os.path.join(base_outd, "TRB", os.path.basename(cols[11]))},
                "TRG": {"src": cols[13], "out": os.path.join(base_outd, "TRG", os.path.basename(cols[13]))},
                "TRD": {"src": cols[15], "out": os.path.join(base_outd, "TRD", os.path.basename(cols[15]))},
                "TRA": {"src": cols[17], "out": os.path.join(base_outd, "TRA", os.path.basename(cols[17]))}
            }

            for _, info in paths.items():
                os.makedirs(os.path.dirname(info["out"]), exist_ok=True)
            
            standard_args = [
                (cols[2], cols[3], paths["IGK"]["out"]), (cols[4], cols[5], paths["IGL"]["out"]),
                (cols[6], cols[7], paths["IGH"]["out"]), (cols[10], cols[11], paths["TRB"]["out"]),
                (cols[12], cols[13], paths["TRG"]["out"]), (cols[16], cols[17], paths["TRA"]["out"]),
                (cols[14], cols[15], paths["TRD"]["out"]),
            ]
            for g, imp, out in standard_args:
                append_pos_import_genes_internal(g, imp, out)
            ighc_append_pos_internal(cols[8], cols[9], paths["IGHC"]["out"])

def parse_mpileup_and_calculate(lines, total_positions):
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    positions_with_10x = 0
    mismatch_list = []
    match_list = []

    for line in lines:
        if not line.strip(): continue
        parts = line.strip().split()
        if len(parts) < 5: continue
        coverage = int(parts[3])
        read_bases = parts[4]
        total_reads += coverage

        match_count = sum(1 for c in read_bases if c in ['.', ','])
        mismatch_count = coverage - match_count
        mismatch_list.append(str(mismatch_count))
        match_list.append(str(match_count))

        if coverage >= 10: positions_with_10x += 1
        
        mismatch_rate = mismatch_count / coverage if coverage > 0 else 0
        match_rate = match_count / coverage if coverage > 0 else 0
        
        if mismatch_rate > 0.2: mismatched_positions += 1
        if match_rate > 0.8: matched_positions += 1

    avg_reads_per_position = total_reads / total_positions if total_positions > 0 else 0
    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0

    return (total_positions, avg_reads_per_position, mismatched_positions, matched_positions,
            ":".join(mismatch_list), ":".join(match_list), percent_accuracy, positions_with_10x)

def get_read_support_vdj3(fofn, scratch):
    gene_types = ["IGK", "IGL", "IGH", "TRB", "TRG", "TRD", "TRA"]

    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for cols in reader:
            if not cols or len(cols) < 17: continue
            sample = cols[0]
            bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
            ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")
            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")

            if not os.path.isfile(bam_file + ".bai"): safe_run(["samtools", "index", bam_file])

            for gene_type in gene_types:
                import_out = os.path.join(base_outd, gene_type, f"{sample}_make_gene_file_imported.csv")
                if not os.path.isfile(import_out): continue
                
                tmp_file = import_out + "_read_support.tmp"

                with open(import_out, 'r', newline='') as f_in, open(tmp_file, 'w', newline='') as f_tmp:
                    csv_in = csv.DictReader(f_in)
                    csv_out = csv.writer(f_tmp)

                    
                    header_row = [
                        "Total_Positions", "Average_Coverage",
                        "Mismatched_Positions",          # Previously ...Coverage_10_Or_Greater
                        "Matched_Positions",             # Previously ...Coverage_10_Or_Greater
                        "Position_Mismatches", "Position_Matches",
                        "Percent_Accuracy", "Positions_With_At_Least_10x_Coverage",
                        "Fully_Spanning_Reads", "Fully_Spanning_Reads_100%_Match"
                    ]
                    csv_out.writerow(header_row)

                    for row_data in csv_in:
                        contig_raw = row_data.get('contig', '')
                        contig = contig_raw.split(',')[0].strip()
                        gene = row_data.get('gene', '')

                        try:
                            start = int(float(row_data['REGION_start']))
                            end = int(float(row_data['REGION_end']))
                        except (ValueError, KeyError):
                            csv_out.writerow([""]*10)
                            continue
                        
                        region = f"{contig}:{start}-{end}"
                        contig_filename = contig.replace('/', '_')
                        tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{start}_{end}.bam")
                        os.makedirs(os.path.dirname(tmp_bam), exist_ok=True)
                        
                        safe_run(["samtools", "view", "-F", "0x100", "-F", "0x800", "-b", bam_file, "-o", tmp_bam, "-U", "/dev/null", region])
                        safe_run(["samtools", "index", tmp_bam])

                        result = safe_run(["samtools", "mpileup", "-f", ref, "-r", region, tmp_bam])
                        output_fields = [""] * 10 
                        if result and result.returncode == 0:
                            mpileup_lines = result.stdout.split('\n')
                            total_positions = end - start + 1
                            stats = parse_mpileup_and_calculate(mpileup_lines, total_positions)
                            
                            # --- SANITY CHECK: Logic Consistency ---
                            validate_vdj_stats(sample, gene, stats)
                            # ---------------------------------------

                            output_fields = list(stats)

                            # --- SANITY CHECK: Sequence ---
                            sequence = get_sequence_from_row(row_data, gene)
                            if not sequence or len(sequence) < 1:
                                print(f"[WARNING] {sample} {gene}: Sequence missing or short (<1bp).")
                            # ------------------------------

                            if sequence:
                                full_span, perf_match = count_matching_reads_internal(tmp_bam, contig, start, end, sequence)
                                output_fields.append(full_span)
                                output_fields.append(perf_match)
                            else:
                                output_fields.append(0)
                                output_fields.append(0)

                        csv_out.writerow(output_fields)
                        try:
                            os.remove(tmp_bam)
                            os.remove(tmp_bam + ".bai")
                        except OSError: pass

                # --- SANITY CHECK: Merge Verification ---
                try:
                    df_orig = pd.read_csv(import_out)
                    df_new = pd.read_csv(tmp_file)
                    
                    if len(df_orig) != len(df_new):
                        print(f"[CRITICAL] {sample} {gene_type}: Merge mismatch! Orig={len(df_orig)}, New={len(df_new)}. Skipping merge.")
                    else:
                        df_combined = pd.concat([df_orig, df_new], axis=1)
                        if 'subject' in df_combined.columns: df_combined['subject'] = sample
                        if 'sample_name' in df_combined.columns: df_combined['sample_name'] = sample
                        df_combined.to_csv(import_out.replace(".csv", "_with_read_support.csv"), index=False)
                except Exception as e:
                    print(f"[WARNING] Failed merging CSVs for {sample}: {e}")
                
                try: os.remove(tmp_file)
                except OSError: pass

def parse_mpileup_ighc(lines, total_positions):
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    mismatch_list = []
    match_list = []
    mism_lt10 = mism_ge10 = mat_lt10 = mat_ge10 = 0

    for line in lines:
        if not line.strip(): continue
        parts = line.strip().split()
        if len(parts) < 5: continue
        coverage = int(parts[3])
        read_bases = parts[4]
        total_reads += coverage
        
        matches = sum(1 for c in read_bases if c in ['.', ','])
        mismatches = coverage - matches
        mismatch_list.append(str(mismatches))
        match_list.append(str(matches)) 

        mismatch_rate = mismatches / coverage if coverage > 0 else 0
        match_rate = matches / coverage if coverage > 0 else 0

        if mismatch_rate > 0.2:
            mismatched_positions += 1
            if coverage < 10: mism_lt10 += 1
            else: mism_ge10 += 1

        if match_rate >= 0.8:
            matched_positions += 1
            if coverage < 10: mat_lt10 += 1
            else: mat_ge10 += 1

    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0
    avg_cov = (total_reads/total_positions) if total_positions > 0 else 0
    pos_10x = (mism_ge10 + mat_ge10)

    return (total_positions, total_reads, avg_cov, mismatched_positions, matched_positions,
            ":".join(mismatch_list), ":".join(match_list), mism_lt10, mism_ge10, mat_lt10, mat_ge10,
            pos_10x, percent_accuracy)

def get_read_support_ighc(fofn, scratch):
    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for cols in reader:
            if not cols or len(cols) < 17: continue
            sample = cols[0]
            bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
            ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")
            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")
            
            if not os.path.isfile(bam_file + ".bai"): safe_run(["samtools", "index", bam_file])
            gene_type = "IGHC"
            import_out = os.path.join(base_outd, gene_type, f"{sample}_make_gene_file_imported.csv")
            if not os.path.isfile(import_out): continue

            tmp_file = import_out + "_read_support.tmp"
            with open(import_out, 'r', newline='') as f_in, open(tmp_file, 'w', newline='') as f_tmp:
                csv_in = csv.DictReader(f_in)
                csv_out = csv.writer(f_tmp)
                
                header_row = [
                    "Total_Positions", "Total_Reads_by_Positions", "Average_Coverage",
                    "Mismatched_Positions", "Matched_Positions",
                    "Position_Mismatches", "Position_Matches",
                    "Mismatched_Positions_Coverage_Less_Than_10",
                    "Mismatched_Positions_Coverage_10_Or_Greater",
                    "Matched_Positions_Coverage_Less_Than_10",
                    "Matched_Positions_Coverage_10_Or_Greater",
                    "Positions_With_At_Least_10x_Coverage", "Percent_Accuracy",
                    "Fully_Spanning_Reads", "Fully_Spanning_Reads_100%_Match",
                    "Allele_reads_100_Match_e1", "Allele_reads_100_Match_e2", "Allele_reads_100_Match_e3",
                    "Allele_reads_100_Match_e4", "Allele_reads_100_Match_e5", "Allele_reads_100_Match_e6",
                    "Allele_reads_100_Match_e7", "Allele_reads_100_Match_e8", "Allele_reads_100_Match_e9",
                    "Allele_reads_fully_spanning_e1", "Allele_reads_fully_spanning_e2", "Allele_reads_fully_spanning_e3",
                    "Allele_reads_fully_spanning_e4", "Allele_reads_fully_spanning_e5", "Allele_reads_fully_spanning_e6",
                    "Allele_reads_fully_spanning_e7", "Allele_reads_fully_spanning_e8", "Allele_reads_fully_spanning_e9"
                ]
                csv_out.writerow(header_row)

                for row_data in csv_in:
                    contig = row_data.get('contig', '').split(',')[0].strip()
                    contig_filename = contig.replace('/', '_')
                    bed_file = os.path.join(base_outd, gene_type, f"{contig_filename}.bed")
                    os.makedirs(os.path.dirname(bed_file), exist_ok=True)

                    min_start, max_end, total_positions = None, None, 0
                    with open(bed_file, 'w') as bed_out:
                        for i in range(1, 10):
                            s_key, e_key = f"C-EXON_{i}_start", f"C-EXON_{i}_end"
                            if s_key in row_data and e_key in row_data:
                                try:
                                    s, e = int(float(row_data[s_key])) - 1, int(float(row_data[e_key]))
                                    if e > 0:
                                        bed_out.write(f"{contig}\t{s}\t{e}\n")
                                        total_positions += (e - s)
                                        if min_start is None or s < min_start: min_start = s
                                        if max_end is None or e > max_end: max_end = e
                                except ValueError: pass
                    
                    if min_start is None:
                        csv_out.writerow([""] * 33)
                        continue

                    region = f"{contig}:{min_start+1}-{max_end}"
                    tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{min_start+1}_{max_end}.bam")
                    try:
                        safe_run(["samtools", "view", "-F", "0x100", "-F", "0x800", "-b", bam_file, "-o", tmp_bam, "-U", "/dev/null", region])
                        safe_run(["samtools", "index", tmp_bam])
                    except Exception: 
                        csv_out.writerow([""] * 33); continue

                    result = safe_run(["samtools", "mpileup", "-f", ref, "-l", bed_file, tmp_bam])
                    output_fields = [""] * 33
                    if result and result.returncode == 0:
                        lines = result.stdout.strip().split('\n')
                        stats = parse_mpileup_ighc(lines, total_positions)
                        output_fields = list(stats)
                        
                        seqs, exons = get_ighc_sequences_from_row(row_data)
                        full, full_all, perf_matches, perf_spans = count_ighc_matches_internal(tmp_bam, contig, min_start+1, max_end, seqs, exons)
                        
                        output_fields.append(full)
                        output_fields.append(full_all)
                        output_fields.extend((perf_matches + [""] * 9)[:9])
                        output_fields.extend((perf_spans + [""] * 9)[:9])

                    csv_out.writerow(output_fields)
                    try: os.remove(tmp_bam); os.remove(tmp_bam + ".bai")
                    except OSError: pass

            try:
                df_orig = pd.read_csv(import_out)
                df_new = pd.read_csv(tmp_file)
                df_combined = pd.concat([df_orig, df_new], axis=1)
                if 'subject' in df_combined.columns: df_combined['subject'] = sample
                if 'sample_name' in df_combined.columns: df_combined['sample_name'] = sample
                df_combined.to_csv(import_out.replace(".csv", "_with_read_support.csv"), index=False)
            except Exception as e: print(f"[WARNING] IGHC Merge Failed: {e}")
            try: os.remove(tmp_file)
            except OSError: pass

def main():
    if len(sys.argv) < 7:
        print("Usage: python3 script.py <fofn> <reffn> <IG_loci> <threads> <scratch> <minimap_option>")
        sys.exit(1)
    
    fofn, reffn, IG_loci, threads, scratch, minimap_option = sys.argv[1:7]
    masked_ref = os.path.join(scratch, "ref_IG_masked.fasta")
    
    run_make_ref_masked(reffn, IG_loci, scratch)
    run_map_ccs_to_pers(fofn, scratch, masked_ref, minimap_option, threads)
    run_append_pos(fofn, scratch)
    get_read_support_vdj3(fofn, scratch)
    get_read_support_ighc(fofn, scratch)
    
    try: os.remove(masked_ref)
    except OSError: pass

if __name__ == "__main__":
    main()