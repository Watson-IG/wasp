#!/usr/bin/env python3
"""
This script reproduces the logic of the provided Bash script in Python.
It does not halt on errors (similar to not using "set -e" in Bash), so it
will continue processing other samples/regions even if a particular command fails.
AWK-based calculations are re-implemented in Python to preserve the same logic.
Refactored to use the 'csv' module for robust file handling.
"""

import sys
import os
import subprocess
import csv

def safe_run(cmd_args):
    """
    Run a subprocess command without halting on error.
    Prints a warning if the command fails, but continues execution.
    """
    try:
        result = subprocess.run(cmd_args, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[WARNING] Command failed: {' '.join(cmd_args)}")
            print(f"stderr:\n{result.stderr}")
        return result
    except Exception as e:
        print(f"[WARNING] Exception while running command: {' '.join(cmd_args)}")
        print(e)
    return None

def run_make_ref_masked(reffn, IG_loci, scratch):
    """
    Masks the given reference with bedtools using the IG_loci regions,
    then indexes the masked reference with samtools.
    """
    masked_ref = os.path.join(scratch, "ref_IG_masked.fasta")

    # bedtools maskfasta
    safe_run([
        "bedtools", "maskfasta",
        "-fi", reffn,
        "-bed", IG_loci,
        "-fo", masked_ref
    ])
    
    # samtools faidx
    safe_run(["samtools", "faidx", masked_ref])

def run_map_ccs_to_pers(fofn, scratch, mask_ref, minimap_option, threads):
    """
    Reads a tab-delimited FOFN (list of samples and inputs).
    For each sample:
      - Creates output directories
      - Extracts contigs from the assembly BAM into a FASTA
      - Concatenates masked ref + personalized contigs
      - Aligns reads.fasta to that personalized reference
      - Sorts and indexes the alignment
    """
    reads_fasta = os.path.join(scratch, "reads.fasta")

    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for cols in reader:
            if not cols: continue 

            # sample, asm_bam, chr2_gene, chr2_import, chr22_gene, chr22_import,
            # igh_gene, igh_import, ighc_gene, ighc_import,
            # trb_gene, trb_import, trg_gene, trg_import,
            # trd_gene, trd_import, tra_gene, tra_import, ccs_bam
            
            if len(cols) <= 18:
                continue

            sample     = cols[0]
            asm_bam    = cols[1]
            ccs_bam    = cols[18]

            print(f"Sample is {sample}")
            print(f"CCS bam is {ccs_bam}")

            outd = os.path.join(scratch, "read_support", sample)
            os.makedirs(os.path.join(outd, "ccs_to_pers"), exist_ok=True)

            safe_run(["samtools", "faidx", reads_fasta])

            contigs_fa = os.path.join(outd, "ccs_to_pers", "contigs.fasta")
            try:
                with subprocess.Popen(
                    ["samtools", "view", "-F", "0x100", "-F", "0x800", asm_bam],
                    stdout=subprocess.PIPE,
                    text=True
                ) as proc, open(contigs_fa, 'w') as contigs_out:
                    for bam_line in proc.stdout:
                        fields = bam_line.strip().split('\t')
                        if len(fields) < 10:
                            continue
                        read_name = fields[0]
                        seq       = fields[9]
                        contigs_out.write(f">{read_name}\n{seq}\n")
            except Exception as e:
                print(f"[WARNING] Failed extracting contigs from {asm_bam} for sample {sample}")
                print(e)
            
            safe_run(["samtools", "faidx", contigs_fa])
            
            pers_ref = os.path.join(outd, "ccs_to_pers", "pers_ref.fasta")
            try:
                with open(pers_ref, 'w') as p_out, \
                     open(mask_ref, 'r') as mask_in, \
                     open(contigs_fa, 'r') as contigs_in:
                    p_out.write(mask_in.read())
                    p_out.write(contigs_in.read())
            except Exception as e:
                print(f"[WARNING] Could not create personalized reference for sample {sample}")
                print(e)

            safe_run(["samtools", "faidx", pers_ref])
            
            sam_out = os.path.join(outd, "ccs_to_pers", "output.sam")
            bam_out = os.path.join(outd, "ccs_to_pers", "output.bam")
            sorted_bam = os.path.join(outd, "ccs_to_pers", "output.sorted.bam")

            safe_run([
                "minimap2",
                "-ax", minimap_option,
                "--secondary=yes",
                "-t", str(threads),
                "-L",
                pers_ref,
                reads_fasta,
                "-o", sam_out
            ])
            safe_run(["samtools", "view", "-Sbh", sam_out, "-o", bam_out])
            safe_run(["samtools", "sort", "-@", str(threads), bam_out, "-o", sorted_bam])
            safe_run(["samtools", "index", sorted_bam])

            try:
                os.remove(sam_out)
            except OSError:
                pass

def run_append_pos(fofn, scratch):
    """
    For each sample in the FOFN, calls external Python scripts to append
    positional information to various gene imports.
    """
    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        
        for cols in reader:
            if not cols: continue
            
            sample       = cols[0]
            if len(cols) < 18:
                continue

            chr2_gene    = cols[2]
            chr2_import  = cols[3]
            chr22_gene   = cols[4]
            chr22_import = cols[5]
            igh_gene     = cols[6]
            igh_import   = cols[7]
            ighc_gene    = cols[8]
            ighc_import  = cols[9]
            trb_gene     = cols[10]
            trb_import   = cols[11]
            trg_gene     = cols[12]
            trg_import   = cols[13]
            trd_gene     = cols[14]
            trd_import   = cols[15]
            tra_gene     = cols[16]
            tra_import   = cols[17]

            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")

            paths = {
                "IGK": {
                    "import_src": chr2_import,
                    "import_out": os.path.join(base_outd, "IGK", os.path.basename(chr2_import))
                },
                "IGL": {
                    "import_src": chr22_import,
                    "import_out": os.path.join(base_outd, "IGL", os.path.basename(chr22_import))
                },
                "IGH": {
                    "import_src": igh_import,
                    "import_out": os.path.join(base_outd, "IGH", os.path.basename(igh_import))
                },
                "IGHC": {
                    "import_src": ighc_import,
                    "import_out": os.path.join(base_outd, "IGHC", os.path.basename(ighc_import))
                },
                "TRB": {
                    "import_src": trb_import,
                    "import_out": os.path.join(base_outd, "TRB", os.path.basename(trb_import))
                },
                "TRG": {
                    "import_src": trg_import,
                    "import_out": os.path.join(base_outd, "TRG", os.path.basename(trg_import))
                },
                "TRD": {
                    "import_src": trd_import,
                    "import_out": os.path.join(base_outd, "TRD", os.path.basename(trd_import))
                },
                "TRA": {
                    "import_src": tra_import,
                    "import_out": os.path.join(base_outd, "TRA", os.path.basename(cols[17]))
                }
            }

            os.makedirs(base_outd, exist_ok=True)
            
            for region_type, info in paths.items():
                os.makedirs(os.path.dirname(info["import_out"]), exist_ok=True)
            
            try:
                script_base = "/opt/wasp/scripts/annotation/read-support"
                
                standard_args = [
                    (chr2_gene, chr2_import, paths["IGK"]["import_out"]),
                    (chr22_gene, chr22_import, paths["IGL"]["import_out"]),
                    (igh_gene, igh_import, paths["IGH"]["import_out"]),
                    (trb_gene, trb_import, paths["TRB"]["import_out"]),
                    (trg_gene, trg_import, paths["TRG"]["import_out"]),
                    (tra_gene, tra_import, paths["TRA"]["import_out"]),
                    (trd_gene, trd_import, paths["TRD"]["import_out"]),
                ]
                
                for g, imp, out in standard_args:
                    safe_run([
                        "/opt/wasp/conda/bin/python",
                        f"{script_base}/append_pos_import_genes.py",
                        g, imp, out
                    ])
                
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    f"{script_base}/ighc_append_pos.py",
                    ighc_gene, ighc_import, paths["IGHC"]["import_out"]
                ])
            except Exception as e:
                print(f"[WARNING] Failed to run append_pos for sample {sample}")
                print(e)

def parse_mpileup_and_calculate(lines, total_positions):
    """
    Reproduce the AWK logic from get_read_support_vdj3.
    """
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    positions_with_10x = 0
    mismatch_list = []
    match_list = []

    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        read_bases = parts[4]
        coverage = int(parts[3])
        total_reads += coverage

        match_count = sum(1 for c in read_bases if c in ['.', ','])
        mismatch_count = coverage - match_count

        mismatch_list.append(str(mismatch_count))
        match_list.append(str(match_count))

        if coverage >= 10:
            positions_with_10x += 1
        
        mismatch_rate = mismatch_count / coverage if coverage > 0 else 0
        match_rate = match_count / coverage if coverage > 0 else 0
        
        if mismatch_rate > 0.2:
            mismatched_positions += 1
        if match_rate > 0.8:
            matched_positions += 1

    avg_reads_per_position = total_reads / total_positions if total_positions > 0 else 0
    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0

    mismatch_str = ":".join(mismatch_list)
    match_str = ":".join(match_list)

    return (
        total_positions,
        avg_reads_per_position,
        mismatched_positions,
        matched_positions,
        mismatch_str,
        match_str,
        percent_accuracy,
        positions_with_10x
    )

def get_read_support_vdj3(fofn, scratch):
    """
    For each sample, processes each gene_type in [IGK, IGL, IGH, TRB, TRG, TRD, TRA],
    runs samtools mpileup, and merges results using csv module.
    """
    gene_types = ["IGK", "IGL", "IGH", "TRB", "TRG", "TRD", "TRA"]

    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        
        for cols in reader:
            if not cols or len(cols) < 17:
                continue
            
            sample = cols[0]
            bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
            ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")
            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")

            if not os.path.isfile(bam_file + ".bai"):
                safe_run(["samtools", "index", bam_file])

            for gene_type in gene_types:
                import_out = os.path.join(base_outd, gene_type, f"{sample}_make_gene_file_imported.csv")
                if not os.path.isfile(import_out):
                    continue
                
                tmp_file = import_out + "_read_support.tmp"

                with open(import_out, 'r', newline='') as f_in, \
                     open(tmp_file, 'w', newline='') as f_tmp:
                    
                    csv_in = csv.reader(f_in)
                    csv_out = csv.writer(f_tmp)

                    try:
                        headers = next(csv_in)
                    except StopIteration:
                        continue 

                    try:
                        contig_col = headers.index("contig")
                        start_col = headers.index("REGION_start")
                        end_col = headers.index("REGION_end")
                        gene_col = headers.index("gene")
                    except ValueError:
                        print(f"[WARNING] Missing required columns in {import_out}")
                        continue
                    
                    header_row = [
                        "Total_Positions", "Average_Coverage",
                        "Mismatched_Positions_Coverage_10_Or_Greater",
                        "Matched_Positions_Coverage_10_Or_Greater",
                        "Position_Mismatches", "Position_Matches",
                        "Percent_Accuracy", "Positions_With_At_Least_10x_Coverage",
                        "Fully_Spanning_Reads", "Fully_Spanning_Reads_100%_Match"
                    ]
                    csv_out.writerow(header_row)

                    for row_data in csv_in:
                        if not row_data: continue
                        
                        # CHANGE: Split contig on comma and take the first one
                        contig = row_data[contig_col].split(',')[0].strip()
                        
                        try:
                            start = int(float(row_data[start_col]))
                            end = int(float(row_data[end_col]))
                        except ValueError:
                            continue
                        
                        region = f"{contig}:{start}-{end}"
                        gene = row_data[gene_col]

                        contig_filename = contig.replace('/', '_')
                        tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{start}_{end}.bam")
                        os.makedirs(os.path.dirname(tmp_bam), exist_ok=True)

                        print(f"Processing region {region} for sample {sample}, gene_type={gene_type}")
                        
                        safe_run([
                            "samtools", "view", "-F", "0x100", "-F", "0x800",
                            "-b", bam_file,
                            "-o", tmp_bam, 
                            "-U", "/dev/null",
                            region
                        ])
                        safe_run(["samtools", "index", tmp_bam])

                        mpileup_cmd = [
                            "samtools", "mpileup",
                            "-f", ref,
                            "-r", region,
                            tmp_bam
                        ]
                        result = safe_run(mpileup_cmd)
                        
                        output_fields = [""] * 10 

                        if result and result.returncode == 0:
                            mpileup_lines = result.stdout.split('\n')
                            total_positions = end - start + 1

                            (
                                tpos, avg_cov, mism_pos, mat_pos,
                                mism_str, mat_str, pct_acc, pos_10x
                            ) = parse_mpileup_and_calculate(mpileup_lines, total_positions)

                            match_cmd = [
                                "/opt/wasp/conda/bin/python",
                                "/opt/wasp/scripts/annotation/read-support/match_subsequences.py",
                                tmp_bam, contig, str(start), str(end), gene, import_out
                            ]
                            match_res = safe_run(match_cmd)
                            match_out = ""
                            if match_res and match_res.returncode == 0:
                                match_out = match_res.stdout.strip()

                            output_fields = [
                                tpos, avg_cov, mism_pos, mat_pos,
                                mism_str, mat_str, pct_acc, pos_10x
                            ]
                            
                            if match_out:
                                output_fields.extend(match_out.split(','))
                        
                        csv_out.writerow(output_fields)

                        try:
                            os.remove(tmp_bam)
                            os.remove(tmp_bam + ".bai")
                        except OSError:
                            pass

                combined_file = import_out.replace(".csv", "_combined.csv")
                final_output = import_out.replace(".csv", "_with_read_support.csv")
                try:
                    with open(import_out, 'r', newline='') as f1, \
                         open(tmp_file, 'r', newline='') as f2, \
                         open(combined_file, 'w', newline='') as fc:
                        
                        r1 = csv.reader(f1)
                        r2 = csv.reader(f2)
                        w = csv.writer(fc)
                        
                        for row1, row2 in zip(r1, r2):
                            w.writerow(row1 + row2)
                            
                except Exception as e:
                    print(f"[WARNING] Failed merging {import_out} and {tmp_file}")
                    print(e)
                    continue

                try:
                    with open(combined_file, 'r', newline='') as f_in, \
                         open(final_output, 'w', newline='') as f_out:
                        
                        r = csv.reader(f_in)
                        w = csv.writer(f_out)
                        
                        first_line = True
                        for row in r:
                            if first_line:
                                w.writerow(row)
                                first_line = False
                            else:
                                if len(row) > 2:
                                    row[1] = sample
                                    row[2] = sample
                                w.writerow(row)
                except Exception as e:
                    print(f"[WARNING] Failed updating subject/sample columns in {combined_file}")
                    print(e)
                    continue

                try:
                    os.remove(combined_file)
                    os.remove(tmp_file)
                except OSError:
                    pass

def parse_mpileup_ighc(lines, total_positions):
    """
    Reproduce the AWK logic from get_read_support_ighc.
    """
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    mismatch_list = []
    match_list = []
    mismatched_positions_cov_lt10 = 0
    mismatched_positions_cov_ge10 = 0
    matched_positions_cov_lt10 = 0
    matched_positions_cov_ge10 = 0

    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 5:
            continue
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
            if coverage < 10:
                mismatched_positions_cov_lt10 += 1
            else:
                mismatched_positions_cov_ge10 += 1

        if match_rate >= 0.8:
            matched_positions += 1
            if coverage < 10:
                matched_positions_cov_lt10 += 1
            else:
                matched_positions_cov_ge10 += 1

    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0
    avg_cov = (total_reads/total_positions) if total_positions > 0 else 0
    pos_10x = (mismatched_positions_cov_ge10 + matched_positions_cov_ge10)

    return (
        total_positions, total_reads, avg_cov, mismatched_positions, matched_positions,
        ":".join(mismatch_list), ":".join(match_list),
        mismatched_positions_cov_lt10, mismatched_positions_cov_ge10,
        matched_positions_cov_lt10, matched_positions_cov_ge10,
        pos_10x, percent_accuracy
    )

def get_read_support_ighc(fofn, scratch):
    """
    Specialized read-support calculation for IGHC genes using csv module.
    """
    with open(fofn, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        
        for cols in reader:
            if not cols or len(cols) < 17:
                continue
            
            sample = cols[0]
            bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
            ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")
            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")

            if not os.path.isfile(bam_file + ".bai"):
                safe_run(["samtools", "index", bam_file])

            gene_type = "IGHC"
            import_out = os.path.join(base_outd, gene_type, f"{sample}_make_gene_file_imported.csv")
            if not os.path.isfile(import_out):
                continue

            tmp_file = import_out + "_read_support.tmp"
            
            with open(import_out, 'r', newline='') as f_in, \
                 open(tmp_file, 'w', newline='') as f_tmp:
                
                csv_in = csv.reader(f_in)
                csv_out = csv.writer(f_tmp)

                try:
                    headers = next(csv_in)
                except StopIteration:
                    continue

                try:
                    contig_col = headers.index("contig")
                    gene_col = headers.index("gene")
                except ValueError:
                    continue

                exon_cols = {}
                for i, h in enumerate(headers):
                    if h.startswith("C-EXON_") and (h.endswith("_start") or h.endswith("_end")):
                        exon_cols[h] = i

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
                    if not row_data: continue
                    
                    # CHANGE: Split contig on comma and take the first one
                    contig = row_data[contig_col].split(',')[0].strip()
                    gene = row_data[gene_col]
                    
                    contig_filename = contig.replace('/', '_')
                    bed_file = os.path.join(base_outd, gene_type, f"{contig_filename}.bed")
                    os.makedirs(os.path.dirname(bed_file), exist_ok=True)

                    with open(bed_file, 'w') as bed_out:
                        for k, idx in exon_cols.items():
                            if "_start" in k:
                                exon_start_str = row_data[idx]
                                exon_end_key = k.replace("_start", "_end")
                                if exon_end_key not in exon_cols:
                                    continue
                                exon_end_str = row_data[exon_cols[exon_end_key]]

                                try:
                                    exon_start = int(float(exon_start_str)) - 1
                                    exon_end = int(float(exon_end_str))
                                except ValueError:
                                    continue
                                
                                if exon_end > 0:
                                    bed_out.write(f"{contig}\t{exon_start}\t{exon_end}\n")
                    
                    total_positions = 0
                    min_start = None
                    max_end = None
                    try:
                        with open(bed_file, 'r') as bed_in:
                            for b_line in bed_in:
                                parts = b_line.strip().split('\t')
                                if len(parts) < 3: continue
                                s = int(parts[1])
                                e = int(parts[2])
                                total_positions += (e - s)
                                
                                if min_start is None or s < min_start: min_start = s
                                if max_end is None or e > max_end: max_end = e
                    except Exception as e:
                        print(f"[WARNING] Error reading bed file {bed_file}: {e}")
                        csv_out.writerow([""] * 33)
                        continue

                    if min_start is None or max_end is None:
                        csv_out.writerow([""] * 33)
                        continue

                    region_start = min_start + 1
                    region_end = max_end
                    region_str = f"{contig}:{region_start}-{region_end}"

                    tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{region_start}_{region_end}.bam")
                    try:
                        safe_run([
                            "samtools", "view",
                            "-F", "0x100", "-F", "0x800",
                            "-b", bam_file,
                            "-o", tmp_bam,
                            "-U", "/dev/null",
                            region_str
                        ])
                        safe_run(["samtools", "index", tmp_bam])
                    except Exception:
                        csv_out.writerow([""] * 33)
                        continue

                    mpileup_cmd = [
                        "samtools", "mpileup",
                        "-f", ref,
                        "-l", bed_file,
                        tmp_bam
                    ]
                    result = safe_run(mpileup_cmd)
                    
                    output_fields = [""] * 33

                    if result and result.returncode == 0:
                        mpileup_lines = result.stdout.strip().split('\n')
                        (
                            tpos, tot_reads, avg_cv, mism_pos, mat_pos,
                            mism_str, mat_str,
                            mism_lt10, mism_ge10, mat_lt10, mat_ge10,
                            pos_10x, pct_acc
                        ) = parse_mpileup_ighc(mpileup_lines, total_positions)

                        match_cmd = [
                            "/opt/wasp/conda/bin/python",
                            "/opt/wasp/scripts/annotation/read-support/ighc_match3.py",
                            tmp_bam, contig, gene, import_out
                        ]
                        match_res = safe_run(match_cmd)
                        match_out = ""
                        if match_res and match_res.returncode == 0:
                            match_out = match_res.stdout.strip()

                        output_fields = [
                            tpos, tot_reads, avg_cv, mism_pos, mat_pos,
                            mism_str, mat_str,
                            mism_lt10, mism_ge10, mat_lt10, mat_ge10,
                            pos_10x, pct_acc
                        ]
                        if match_out:
                            output_fields.extend(match_out.split(','))

                    csv_out.writerow(output_fields)

                    try:
                        os.remove(tmp_bam)
                        os.remove(tmp_bam + ".bai")
                    except OSError:
                        pass

            combined_file = import_out.replace(".csv", "_combined.csv")
            final_output = import_out.replace(".csv", "_with_read_support.csv")
            try:
                with open(import_out, 'r', newline='') as f1, \
                     open(tmp_file, 'r', newline='') as f2, \
                     open(combined_file, 'w', newline='') as fc:
                    
                    r1 = csv.reader(f1)
                    r2 = csv.reader(f2)
                    w = csv.writer(fc)
                    
                    for row1, row2 in zip(r1, r2):
                        w.writerow(row1 + row2)
            except Exception as e:
                print(f"[WARNING] Failed merging IGHC {import_out} and {tmp_file}")
                print(e)
                continue

            try:
                with open(combined_file, 'r', newline='') as f_in, \
                     open(final_output, 'w', newline='') as f_out:
                    
                    r = csv.reader(f_in)
                    w = csv.writer(f_out)
                    
                    first_line = True
                    for row in r:
                        if first_line:
                            w.writerow(row)
                            first_line = False
                        else:
                            if len(row) > 2:
                                row[1] = sample
                                row[2] = sample
                            w.writerow(row)
            except Exception as e:
                print(f"[WARNING] Failed updating subject/sample columns in {combined_file}")
                print(e)
                continue

            try:
                os.remove(combined_file)
                os.remove(tmp_file)
            except OSError:
                pass

def main():
    if len(sys.argv) < 7:
        print("Usage: python3 this_script.py <fofn> <reffn> <IG_loci> <threads> <scratch> <minimap_option>")
        sys.exit(1)

    fofn = sys.argv[1]
    reffn = sys.argv[2]
    IG_loci = sys.argv[3]
    threads = sys.argv[4]
    scratch = sys.argv[5]
    minimap_option = sys.argv[6]

    masked_ref = os.path.join(scratch, "ref_IG_masked.fasta")

    run_make_ref_masked(reffn, IG_loci, scratch)
    run_map_ccs_to_pers(fofn, scratch, masked_ref, minimap_option, threads)
    run_append_pos(fofn, scratch)
    get_read_support_vdj3(fofn, scratch)
    get_read_support_ighc(fofn, scratch)

    try:
        os.remove(masked_ref)
    except OSError:
        pass

if __name__ == "__main__":
    main()