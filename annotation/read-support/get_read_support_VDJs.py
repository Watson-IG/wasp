#!/usr/bin/env python3
"""
This script reproduces the logic of the provided Bash script in Python.
It does not halt on errors (similar to not using "set -e" in Bash), so it
will continue processing other samples/regions even if a particular command fails.
AWK-based calculations are re-implemented in Python to preserve the same logic.
"""

import sys
import os
import subprocess

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
    # Path to the HiFi reads in FASTA form (assumes this is a common file per the original script logic)
    reads_fasta = os.path.join(scratch, "reads.fasta")

    with open(fofn, 'r') as infile:
        for line in infile:
            # sample, asm_bam, chr2_gene, chr2_import, chr22_gene, chr22_import,
            # igh_gene, igh_import, ighc_gene, ighc_import,
            # trb_gene, trb_import, trg_gene, trg_import,
            # trd_gene, trd_import, tra_gene, tra_import, ccs_bam  (
            cols = line.strip().split('\t')
            
            sample     = cols[0]
            asm_bam    = cols[1]
            ccs_bam    = cols[18]  # The original script references ccs_bam at the end

            print(f"Sample is {sample}")
            print(f"CCS bam is {ccs_bam}")

            outd = os.path.join(scratch, "read_support", sample)
            os.makedirs(os.path.join(outd, "ccs_to_pers"), exist_ok=True)

            # Make sure reads_fasta is indexed
            # The original script calls: samtools faidx ${scratch}/reads.fasta
            safe_run(["samtools", "faidx", reads_fasta])

            # Extract contigs from the assembly BAM into a FASTA
            # Equivalent to: samtools view -F 0x100 -F 0x800 asm_bam | awk '{print ">"$1"\n"$10}'
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
            
            # Index the new contigs.fasta
            safe_run(["samtools", "faidx", contigs_fa])
            
            # Create the personalized reference
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

            # Index personalized ref
            safe_run(["samtools", "faidx", pers_ref])
            
            # Align reads.fasta to personalized ref
            sam_out = os.path.join(outd, "ccs_to_pers", "output.sam")
            bam_out = os.path.join(outd, "ccs_to_pers", "output.bam")
            sorted_bam = os.path.join(outd, "ccs_to_pers", "output.sorted.bam")

            # Minimap2 alignment
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
            # Convert SAM to BAM
            safe_run(["samtools", "view", "-Sbh", sam_out, "-o", bam_out])
            # Sort
            safe_run(["samtools", "sort", "-@", str(threads), bam_out, "-o", sorted_bam])
            # Index
            safe_run(["samtools", "index", sorted_bam])

            # Clean up
            try:
                os.remove(sam_out)
            except OSError:
                pass

def run_append_pos(fofn, scratch):
    """
    For each sample in the FOFN, calls external Python scripts to append
    positional information to various gene imports.
    """
    with open(fofn, 'r') as infile:
        for line in infile:
            cols = line.strip().split()
            
            sample       = cols[0]
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
            # The next columns could be ccs_bam or so, but we only need these for appending.

            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")

            # Prepare subdirectories and outputs
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
                    "import_src": tra_import,  # Actually the script uses tra_gene, tra_import
                    "import_out": os.path.join(base_outd, "TRA", os.path.basename(cols[17]) if len(cols) > 17 else "tra_import.csv")
                }
            }

        

            # Make directories and run external append scripts
            os.makedirs(base_outd, exist_ok=True)
            
            # The original script calls these python commands for the first seven gene sets
            # and a separate script for IGHC. We'll replicate that:
            # /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py
            for region_type, info in paths.items():
                os.makedirs(os.path.dirname(info["import_out"]), exist_ok=True)
            
            try:
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    chr2_gene, chr2_import, paths["IGK"]["import_out"]
                ])
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    chr22_gene, chr22_import, paths["IGL"]["import_out"]
                ])
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    igh_gene, igh_import, paths["IGH"]["import_out"]
                ])
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    trb_gene, trb_import, paths["TRB"]["import_out"]
                ])
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    trg_gene, trg_import, paths["TRG"]["import_out"]
                ])
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    tra_gene, tra_import, paths["TRA"]["import_out"]
                ])
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py",
                    trd_gene, trd_import, paths["TRD"]["import_out"]
                ])
                
                # IGHC uses a different script
                safe_run([
                    "/opt/wasp/conda/bin/python",
                    "/opt/wasp/scripts/annotation/read-support/ighc_append_pos.py",
                    ighc_gene, ighc_import, paths["IGHC"]["import_out"]
                ])
            except Exception as e:
                print(f"[WARNING] Failed to run append_pos for sample {sample}")
                print(e)

def parse_mpileup_and_calculate(lines, total_positions):
    """
    Reproduce the AWK logic from get_read_support_vdj3:
      - total_reads = sum of coverage
      - coverage = len(column5)
      - mismatches = number of non-. or , in column5
      - matches = coverage - mismatches
      - mismatch_rate > 0.2 => mismatched_positions++
      - match_rate > 0.8 => matched_positions++
      - coverage >= 10 => positions_with_10x++
      - mismatch_list: colon-separated mismatch counts
      - match_list: colon-separated match counts
      - avg_reads_per_position = total_reads / total_positions
      - percent_accuracy = (matched_positions / total_positions) * 100
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
        coverage = len(read_bases)
        total_reads += coverage

        # Count mismatches
        # AWK used gensub(/[.,]/, "", "g", $5) => remove '.' and ',' => everything left is mismatches
        mismatch_count = sum(1 for c in read_bases if c not in ['.', ','])
        match_count = coverage - mismatch_count  # or sum(1 for c in read_bases if c in ['.', ','])

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

    # Format the mismatch_list and match_list with colons
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
    runs samtools mpileup, re-implements the AWK logic in Python,
    and merges the result with the external match_subsequences.py output.
    """
    gene_types = ["IGK", "IGL", "IGH", "TRB", "TRG", "TRD", "TRA"]

    with open(fofn, 'r') as infile:
        for line in infile:
            cols = line.strip().split()
            if len(cols) < 17:
                continue
            
            sample = cols[0]
            bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
            ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")
            base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")

            # Index if needed
            if not os.path.isfile(bam_file + ".bai"):
                safe_run(["samtools", "index", bam_file])

            for gene_type in gene_types:
                import_out = os.path.join(base_outd, gene_type, f"{sample}_make_gene_file_imported.csv")
                if not os.path.isfile(import_out):
                    continue
                
                tmp_file = import_out + "_read_support.tmp"
                # Write a header to the tmp_file (to match the script's initial CSV header)
                with open(tmp_file, 'w') as f:
                    f.write("Total_Positions,Average_Coverage,Mismatched_Positions_Coverage_10_Or_Greater,Matched_Positions_Coverage_10_Or_Greater,") #! match, mismatch are not checked for 10x
                    f.write("Position_Mismatches,Position_Matches,Percent_Accuracy,")
                    f.write("Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,")
                    f.write("Fully_Spanning_Reads_100%_Match\n")

                # Identify column indices for contig, REGION_start, REGION_end, gene
                with open(import_out, 'r') as f_in:
                    header_line = f_in.readline().strip()
                    headers = header_line.split(',')
                    contig_col = None
                    start_col = None
                    end_col = None
                    gene_col = None

                    for i, h in enumerate(headers):
                        if h == "contig":
                            contig_col = i
                        elif h == "REGION_start":
                            start_col = i
                        elif h == "REGION_end":
                            end_col = i
                        elif h == "gene":
                            gene_col = i
                
                # Process each line after the header
                with open(import_out, 'r') as f_in, open(tmp_file, 'a') as f_tmp:
                    next(f_in)  # skip header
                    for row in f_in:
                        row_data = row.strip().split(',')
                        if len(row_data) <= max(contig_col, start_col, end_col, gene_col):
                            continue
                        
                        contig = row_data[contig_col]
                        try:
                            start = int(float(row_data[start_col]))
                            end = int(float(row_data[end_col]))
                        except ValueError:
                            continue
                        region = f"{contig}:{start}-{end}"
                        gene = row_data[gene_col]

                        # Create temporary region-specific BAM
                        contig_filename = contig.replace('/', '_')
                        tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{start}_{end}.bam")
                        os.makedirs(os.path.dirname(tmp_bam), exist_ok=True)

                        print(f"Processing region {region} for sample {sample}, gene_type={gene_type}")
                        
                        # samtools view -F 0x100 -F 0x800 -b ...
                        safe_run([
                            "samtools", "view", "-F", "0x100", "-F", "0x800",
                            "-b", bam_file,
                            "-o", tmp_bam, 
                            "-U", "/dev/null",
                            region
                        ])
                        safe_run(["samtools", "index", tmp_bam])

                        # Run mpileup
                        mpileup_cmd = [
                            "samtools", "mpileup",
                            "-f", ref,
                            "-r", region,
                            tmp_bam
                        ]
                        result = safe_run(mpileup_cmd)
                        if not result or result.returncode != 0:
                            # Could not get mpileup output
                            # Write placeholder row to keep alignment with final output
                            f_tmp.write(",,,,,,,,,\n")
                            try:
                                os.remove(tmp_bam)
                                os.remove(tmp_bam + ".bai")
                            except OSError:
                                pass
                            continue
                        
                        mpileup_lines = result.stdout.split('\n')
                        total_positions = end - start + 1

                        # Reproduce the AWK logic in Python
                        (
                            tpos,
                            avg_cov,
                            mismatched_positions,
                            matched_positions,
                            mismatch_str,
                            match_str,
                            pct_acc,
                            pos_10x
                        ) = parse_mpileup_and_calculate(mpileup_lines, total_positions)

                        # Now call the external match_subsequences.py to capture its stdout
                        match_cmd = [
                            "/opt/wasp/conda/bin/python",
                            "/opt/wasp/scripts/annotation/read-support/match_subsequences.py",
                            tmp_bam,
                            contig,
                            str(start),
                            str(end),
                            gene,
                            import_out
                        ]
                        match_res = safe_run(match_cmd)
                        match_out = ""
                        if match_res and match_res.returncode == 0:
                            match_out = match_res.stdout.strip()
                        
                        # Combine and write to tmp_file
                        row_str = f"{tpos},{avg_cov},{mismatched_positions},{matched_positions},{mismatch_str},{match_str},{pct_acc},{pos_10x}"
                        if match_out:
                            row_str = row_str + "," + match_out
                        f_tmp.write(row_str + "\n")

                        # Cleanup region-specific BAM
                        try:
                            os.remove(tmp_bam)
                            os.remove(tmp_bam + ".bai")
                        except OSError:
                            pass

                # Merge the newly created tmp_file with the original import_out
                combined_file = import_out.replace(".csv", "_combined.csv")
                final_output = import_out.replace(".csv", "_with_read_support.csv")
                try:
                    with open(import_out, 'r') as f1, open(tmp_file, 'r') as f2, open(combined_file, 'w') as fc:
                        for l1, l2 in zip(f1, f2):
                            fc.write(l1.strip() + "," + l2.strip() + "\n")
                except Exception as e:
                    print(f"[WARNING] Failed merging {import_out} and {tmp_file}")
                    print(e)
                    continue

                # Update columns 2 and 3 (Subject, Sample_Name) with the sample name
                try:
                    with open(combined_file, 'r') as f_in, open(final_output, 'w') as f_out:
                        first_line = True
                        for line_comb in f_in:
                            line_comb = line_comb.strip()
                            if first_line:
                                f_out.write(line_comb + "\n")
                                first_line = False
                            else:
                                parts = line_comb.split(',')
                                if len(parts) > 2:
                                    parts[1] = sample
                                    parts[2] = sample
                                f_out.write(",".join(parts) + "\n")
                except Exception as e:
                    print(f"[WARNING] Failed updating subject/sample columns in {combined_file}")
                    print(e)
                    continue

                # Cleanup
                try:
                    os.remove(combined_file)
                    os.remove(tmp_file)
                except OSError:
                    pass

def parse_mpileup_ighc(lines, total_positions):
    """
    Reproduce the AWK logic from get_read_support_ighc:
      - total_reads (sum of coverage)
      - mismatched_positions
      - matched_positions
      - mismatch_list, match_list (colon-separated)
      - coverage-based subcounts:
          mismatched_positions_coverage_less_than_10
          mismatched_positions_coverage_10_or_greater
          matched_positions_coverage_less_than_10
          matched_positions_coverage_10_or_greater
      - percent_accuracy = (matched_positions/total_positions) * 100
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
    avg_cov = (total_reads/total_positions)
    pos_10x = (mismatched_positions_cov_ge10 + matched_positions_cov_ge10)

    return (
        total_positions,
        total_reads,
        avg_cov,
        mismatched_positions,
        matched_positions,
        ":".join(mismatch_list),
        ":".join(match_list),
        mismatched_positions_cov_lt10,
        mismatched_positions_cov_ge10,
        matched_positions_cov_lt10,
        matched_positions_cov_ge10,
        pos_10x,
        percent_accuracy
    )

def get_read_support_ighc(fofn, scratch):
    """
    Specialized read-support calculation for IGHC genes.
    Splits up exons, merges coverage stats with an external script,
    and appends results to the existing CSV.
    """
    with open(fofn, 'r') as infile:
        for line in infile:
            cols = line.strip().split()
            if len(cols) < 17:
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
            
            # Write header to tmp_file
            with open(tmp_file, 'w') as f:
                f.write("Total_Positions,Total_Reads_by_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,")
                f.write("Position_Mismatches,Position_Matches,")
                f.write("Mismatched_Positions_Coverage_Less_Than_10,")
                f.write("Mismatched_Positions_Coverage_10_Or_Greater,")
                f.write("Matched_Positions_Coverage_Less_Than_10,")
                f.write("Matched_Positions_Coverage_10_Or_Greater,")
                f.write("Positions_With_At_Least_10x_Coverage,Percent_Accuracy,")
                f.write("Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match,")
                f.write("Allele_reads_100_Match_e1,Allele_reads_100_Match_e2,Allele_reads_100_Match_e3,")
                f.write("Allele_reads_100_Match_e4,Allele_reads_100_Match_e5,Allele_reads_100_Match_e6,")
                f.write("Allele_reads_100_Match_e7,Allele_reads_100_Match_e8,Allele_reads_100_Match_e9,")
                f.write("Allele_reads_fully_spanning_e1,Allele_reads_fully_spanning_e2,Allele_reads_fully_spanning_e3,")
                f.write("Allele_reads_fully_spanning_e4,Allele_reads_fully_spanning_e5,Allele_reads_fully_spanning_e6,")
                f.write("Allele_reads_fully_spanning_e7,Allele_reads_fully_spanning_e8,Allele_reads_fully_spanning_e9\n")

            # Identify columns we need: contig, gene, plus C-EXON_X_start/end
            with open(import_out, 'r') as f_in:
                header_line = f_in.readline().strip()
                headers = header_line.split(',')
                contig_col = None
                gene_col = None
                exon_cols = {}

                for i, h in enumerate(headers):
                    if h == "contig":
                        contig_col = i
                    elif h == "gene":
                        gene_col = i
                    # We'll store exon columns in a dict
                    # e.g. h might be C-EXON_1_start or C-EXON_1_end
                    if h.startswith("C-EXON_") and (h.endswith("_start") or h.endswith("_end")):
                        exon_cols[h] = i

            # Process each line, build a .bed of all exons, then mpileup on them
            with open(import_out, 'r') as f_in, open(tmp_file, 'a') as f_tmp:
                next(f_in)  # skip header
                for row in f_in:
                    row_data = row.strip().split(',')
                    if len(row_data) <= max(list(exon_cols.values()) + [contig_col, gene_col]):
                        continue
                    
                    contig = row_data[contig_col]
                    gene = row_data[gene_col]
                    
                    # Make a BED file containing all exons for this IGHC entry
                    contig_filename = contig.replace('/', '_')
                    bed_file = os.path.join(base_outd, gene_type, f"{contig_filename}.bed")
                    os.makedirs(os.path.dirname(bed_file), exist_ok=True)

                    # Clear the BED file
                    with open(bed_file, 'w') as bed_out:
                        pass

                    # Write exon lines to the BED
                    # For each EXON_X_start, find the matching EXON_X_end, write contig, start, end
                    for k, idx in exon_cols.items():
                        # k looks like "C-EXON_1_start" or "C-EXON_1_end"
                        if "_start" in k:
                            exon_start_str = row_data[idx]
                            # replace _start with _end to find matching col
                            exon_end_key = k.replace("_start", "_end")
                            if exon_end_key not in exon_cols:
                                continue
                            exon_end_str = row_data[exon_cols[exon_end_key]]

                            try:
                                exon_start = int(float(exon_start_str)) - 1  # zero-based for BED
                                exon_end = int(float(exon_end_str))
                            except ValueError:
                                continue
                            
                            # Only write if exon_end != 0, or valid range
                            if exon_end > 0:
                                # Append to bed file
                                with open(bed_file, 'a') as bed_out:
                                    bed_out.write(f"{contig}\t{exon_start}\t{exon_end}\n")
                    
                    # Now compute the total positions from the BED
                    total_positions = 0
                    try:
                        with open(bed_file, 'r') as bed_in:
                            for b_line in bed_in:
                                parts = b_line.strip().split('\t')
                                if len(parts) < 3:
                                    continue
                                start = int(parts[1])
                                end = int(parts[2])
                                total_positions += (end - start)
                    except Exception as e:
                        print(f"[WARNING] Error reading bed file {bed_file} for sample {sample}: {e}")
                        f_tmp.write("," * 21 + "\n")  # write blank line
                        continue

                  
                    min_start = None
                    max_end = None
                    try:
                        with open(bed_file, 'r') as bed_in:
                            for b_line in bed_in:
                                parts = b_line.strip().split('\t')
                                s = int(parts[1])
                                e = int(parts[2])
                                if min_start is None or s < min_start:
                                    min_start = s
                                if max_end is None or e > max_end:
                                    max_end = e
                    except Exception as e:
                        print(f"[WARNING] Failed scanning bed for sample {sample}")
                        f_tmp.write("," * 21 + "\n")
                        continue

                    if min_start is None or max_end is None:
                        # no exons
                        f_tmp.write("," * 21 + "\n")
                        continue

                    region_start = min_start + 1  # convert back to 1-based for samtools
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
                    except Exception as e:
                        print(f"[WARNING] samtools view/index error for IGHC region {region_str} in sample {sample}")
                        f_tmp.write("," * 21 + "\n")
                        continue

                    # Now run mpileup on the entire bed
                    mpileup_cmd = [
                        "samtools", "mpileup",
                        "-f", ref,
                        "-l", bed_file,
                        tmp_bam
                    ]
                    result = safe_run(mpileup_cmd)
                    if not result or result.returncode != 0:
                        # Could not get mpileup output
                        f_tmp.write("," * 21 + "\n")
                        try:
                            os.remove(tmp_bam)
                            os.remove(tmp_bam + ".bai")
                        except OSError:
                            pass
                        continue

                    mpileup_lines = result.stdout.strip().split('\n')

                    (
                        tpos,
                        total_reads,
                        avg_cov,
                        mismatched_positions,
                        matched_positions,
                        mismatch_str,
                        match_str,
                        mismatched_positions_cov_lt10,
                        mismatched_positions_cov_ge10,
                        matched_positions_cov_lt10,
                        matched_positions_cov_ge10,
                        pos_10x,
                        pct_acc
                    ) = parse_mpileup_ighc(mpileup_lines, total_positions)

                    # External IGHC match
                    match_cmd = [
                        "/opt/wasp/conda/bin/python",
                        "/opt/wasp/scripts/annotation/read-support/ighc_match3.py",
                        tmp_bam,
                        contig,
                        gene,
                        import_out
                    ]
                    match_res = safe_run(match_cmd)
                    match_out = ""
                    if match_res and match_res.returncode == 0:
                        match_out = match_res.stdout.strip()

                    # Write final line
                    row_str = (f"{tpos},{total_reads},{avg_cov},{mismatched_positions},{matched_positions},"
                               f"{mismatch_str},{match_str},"
                               f"{mismatched_positions_cov_lt10},{mismatched_positions_cov_ge10},"
                               f"{matched_positions_cov_lt10},{matched_positions_cov_ge10},"
                               f"{pos_10x},{pct_acc}")
                    if match_out:
                        row_str += "," + match_out
                    f_tmp.write(row_str + "\n")

                    # Cleanup
                    try:
                        os.remove(tmp_bam)
                        os.remove(tmp_bam + ".bai")
                    except OSError:
                        pass

            # Merge the newly created tmp_file with the original import_out
            combined_file = import_out.replace(".csv", "_combined.csv")
            final_output = import_out.replace(".csv", "_with_read_support.csv")
            try:
                with open(import_out, 'r') as f1, open(tmp_file, 'r') as f2, open(combined_file, 'w') as fc:
                    for l1, l2 in zip(f1, f2):
                        fc.write(l1.strip() + "," + l2.strip() + "\n")
            except Exception as e:
                print(f"[WARNING] Failed merging IGHC {import_out} and {tmp_file}")
                print(e)
                continue

            # Update columns 2 and 3 with the sample name
            try:
                with open(combined_file, 'r') as f_in, open(final_output, 'w') as f_out:
                    first_line = True
                    for line_comb in f_in:
                        line_comb = line_comb.strip()
                        if first_line:
                            f_out.write(line_comb + "\n")
                            first_line = False
                        else:
                            parts = line_comb.split(',')
                            if len(parts) > 2:
                                parts[1] = sample
                                parts[2] = sample
                            f_out.write(",".join(parts) + "\n")
            except Exception as e:
                print(f"[WARNING] Failed updating subject/sample columns in {combined_file}")
                print(e)
                continue

            # Cleanup
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

    # Path where we'll store the masked reference
    masked_ref = os.path.join(scratch, "ref_IG_masked.fasta")

    run_make_ref_masked(reffn, IG_loci, scratch)
    run_map_ccs_to_pers(fofn, scratch, masked_ref, minimap_option, threads)
    run_append_pos(fofn, scratch)
    get_read_support_vdj3(fofn, scratch)
    get_read_support_ighc(fofn, scratch)

    # Remove the masked reference
    try:
        os.remove(masked_ref)
    except OSError:
        pass

if __name__ == "__main__":
    main()
