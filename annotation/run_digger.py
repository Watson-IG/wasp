#!/usr/bin/env python3
"""
run_digger.py

Bins assembly contigs by IG/TR locus using BLAST against reference allele
sequences, then runs the digger tool on each per-locus FASTA, and adds
read support to the digger annotation tables using wasptk.

Workflow:
  1. Scan allele_ref_dir for reference FASTAs matching the species prefix
     to discover which loci exist (e.g. IGH, IGK, IGL, TRA, TRB, ...).
  2. Build a combined BLAST database from all non-gapped reference FASTAs.
  3. BLAST the assembly contigs (full_asm_for_digger.fasta) against the DB.
  4. Assign each contig to its best-hit locus.
  5. Write per-locus FASTA files (e.g. IGH_contigs.fasta).
  6. Run digger on each per-locus FASTA, passing the appropriate
     species-specific V/D/J/C reference files.
     - IGK is special: run digger in both + and - sense, then merge outputs.
  7. Map CCS reads to per-locus contigs.
  8. Run wasptk readsupport on each digger output CSV.
"""

import argparse
import glob
import os
import pandas as pd
import re
import subprocess
import sys
from collections import defaultdict, Counter


def discover_loci(allele_ref_dir: str, species: str) -> dict[str, dict[str, str]]:
    """Scan allele_ref_dir and return a dict of locus -> {segment_type: filepath}.

    For example:
        {
            'IGH': {'V': '...IGHV.fasta', 'D': '...IGHD.fasta', ...},
            'IGK': {'V': '...IGKV.fasta', 'J': '...IGKJ.fasta', ...},
            ...
        }

    Only non-gapped FASTAs are included (gapped files are tracked separately).
    """
    # Normalise the species prefix to match filenames (e.g. "human" -> "Homo_sapiens")
    # The reference files use the binomial name, so the caller must pass it correctly.
    pattern = os.path.join(allele_ref_dir, f"{species}_*.fasta")
    files = glob.glob(pattern)

    if not files:
        sys.exit(
            f"ERROR: No reference files found matching '{pattern}'. "
            f"Check -species and -allele_ref_dir values."
        )

    loci: dict[str, dict[str, str]] = defaultdict(dict)
    # Regex to parse e.g. "Homo_sapiens_IGHV.fasta" or "Homo_sapiens_IGHV_gapped.fasta"
    fname_re = re.compile(
        rf"^{re.escape(species)}_([A-Z]{{2,3}})([VDJC])(_gapped)?\.fasta$"
    )

    for fpath in files:
        fname = os.path.basename(fpath)
        m = fname_re.match(fname)
        if not m:
            continue
        locus = m.group(1)       # e.g. "IGH", "IGK", "TRA"
        segment = m.group(2)     # e.g. "V", "D", "J", "C"
        is_gapped = m.group(3)   # e.g. "_gapped" or None

        if is_gapped:
            loci[locus][f"{segment}_gapped"] = fpath
        else:
            loci[locus][segment] = fpath

    if not loci:
        sys.exit(
            f"ERROR: Reference files found but none matched the expected naming "
            f"convention ({species}_<LOCUS><SEGMENT>.fasta)."
        )

    return dict(loci)


def build_blast_db(loci_to_bin: dict[str, dict[str, str]], outdir: str) -> str:
    """Concatenate all non-gapped reference FASTAs and build a BLAST DB.

    Returns the path to the BLAST database.
    """
    db_dir = os.path.join(outdir, "digger_blast_db")
    os.makedirs(db_dir, exist_ok=True)

    combined_fasta = os.path.join(db_dir, "combined_refs.fasta")

    with open(combined_fasta, "w") as out_fh:
        for locus, refs in loci_to_bin.items():
            for segment, ref_file in refs.items():
                if "_gapped" not in segment:
                    with open(ref_file) as in_fh:
                        for line in in_fh:
                            out_fh.write(line)

    db_path = os.path.join(db_dir, "combined_refs")
    subprocess.run(
        ["makeblastdb", "-in", combined_fasta, "-dbtype", "nucl", "-out", db_path],
        check=True,
    )
    return db_path


def blast_contigs(assembly_fasta: str, db_path: str, outdir: str) -> str:
    """BLAST contigs against the combined reference DB.

    Returns path to the BLAST output (outfmt 6).
    """
    blast_out = os.path.join(outdir, "digger_blast_db", "blast_results.tsv")
    subprocess.run(
        [
            "blastn",
            "-query", assembly_fasta,
            "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore slen",
            "-evalue", "1e-5",
            "-max_target_seqs", "20",
            "-out", blast_out,
        ],
        check=True,
    )
    return blast_out


def bin_contigs_by_locus(
    blast_results: str, assembly_fasta: str, loci_to_bin: dict, outdir: str
) -> dict[str, str]:
    """Parse BLAST results, assign each contig to a locus, write per-locus FASTAs.

    TRA and TRD are co-located (TRD is nested within TRA on chr14), so contigs
    hitting either locus are assigned to both.

    Returns a dict of locus -> path to the per-locus contig FASTA.
    """
    # Loci that share the same genomic region — contigs hitting one should
    # also be binned into the other so digger can annotate genes from both.
    COLOCATED_LOCI: dict[str, set[str]] = {
        "TRA": {"TRD"},
        "TRD": {"TRA"},
    }

    # Map each contig to its best-hit locus based on the subject ID.
    # Subject IDs come from the reference FASTAs and contain the locus name,
    # e.g. a hit to "IGHV1-2*01" means the contig belongs to IGH.
    known_loci = set(loci_to_bin.keys())

    # Map each contig to a list of locus hits
    contig_hits: dict[str, list[str]] = defaultdict(list)
    with open(blast_results) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue
            
            qseqid, sseqid = fields[0], fields[1]
            length = int(fields[3])
            slen = int(fields[6])
            
            if length < 0.95 * slen:
                continue

            # Try to extract the locus prefix from the subject sequence ID.
            # Reference sequence names typically start with the locus,
            # e.g. "IGHV1-2*01", "IGKJ5*01", "TRAV1*01"
            for locus in known_loci:
                if locus in sseqid:
                    contig_hits[qseqid].append(locus)
                    break

    contig_to_loci: dict[str, set[str]] = defaultdict(set)
    for qseqid, hits in contig_hits.items():
        if not hits:
            continue
        
        hit_counts = Counter(hits)
        total_hits = len(hits)
        top_locus, top_count = hit_counts.most_common(1)[0]
        
        # If one locus constitutes more than 80% of the valid hits, assign only to that locus.
        # Otherwise, assign to all loci that were hit.
        if (top_count / total_hits) > 0.80:
            assigned_loci = {top_locus}
        else:
            assigned_loci = set(hits)
            
        for locus in assigned_loci:
            contig_to_loci[qseqid].add(locus)
            # Also assign to co-located loci (e.g. TRA <-> TRD)
            for partner in COLOCATED_LOCI.get(locus, set()):
                if partner in known_loci:
                    contig_to_loci[qseqid].add(partner)

    # Read all contigs from the assembly FASTA
    contigs: dict[str, str] = {}
    current_id = None
    current_seq: list[str] = []
    with open(assembly_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    contigs[current_id] = "\n".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = [line.rstrip()]
            else:
                current_seq.append(line.rstrip())
        if current_id is not None:
            contigs[current_id] = "\n".join(current_seq)

    # Write per-locus FASTA files
    binned_dir = os.path.join(outdir, "digger_binned")
    os.makedirs(binned_dir, exist_ok=True)

    locus_fastas: dict[str, str] = {}
    locus_contigs: dict[str, list[str]] = defaultdict(list)

    for contig_id, loci_set in contig_to_loci.items():
        if contig_id in contigs:
            for locus in loci_set:
                locus_contigs[locus].append(contigs[contig_id])

    for locus, seqs in locus_contigs.items():
        fasta_path = os.path.join(binned_dir, f"{locus}_contigs.fasta")
        with open(fasta_path, "w") as fh:
            fh.write("\n".join(seqs) + "\n")
        locus_fastas[locus] = fasta_path

    return locus_fastas


def _build_digger_cmd(
    fasta_path: str,
    output_file: str,
    refs: dict[str, str],
    species: str,
    locus: str,
    sense: str | None = None,
    motif_dir: str | None = None,
) -> list[str]:
    """Build a digger command list."""
    # Map binomial species names used for references to digger's internal motif folder names
    digger_species = species.lower()
    if digger_species == "homo_sapiens":
        digger_species = "human"
    elif digger_species == "macaca_mulatta":
        digger_species = "rhesus_macaque"

    cmd = [
        "digger",
        fasta_path,
        output_file,
        "-species", digger_species,
    ]
    if "V" in refs:
        cmd.extend(["-v_ref", refs["V"]])
    if "D" in refs:
        cmd.extend(["-d_ref", refs["D"]])
    if "J" in refs:
        cmd.extend(["-j_ref", refs["J"]])
    if "V_gapped" in refs:
        cmd.extend(["-v_ref_gapped", refs["V_gapped"]])
    cmd.extend(["-locus", locus])
    if sense is not None:
        cmd.extend(["-sense", sense])
    if motif_dir:
        cmd.extend(["-motif_dir", motif_dir])
    return cmd


def _check_digger_output(csv_path: str, locus: str, sense: str | None = None) -> bool:
    """Check if a digger output CSV has any data rows beyond the header.

    Digger can exit 0 but produce a headers-only file. This catches that.
    Returns True if the table has data rows, False otherwise.
    """
    sense_label = f" (sense {sense})" if sense else ""
    if not os.path.isfile(csv_path):
        print(f"DIGGER_WARNING: {locus}{sense_label} output file not found: {csv_path}")
        return False
    try:
        df = pd.read_csv(csv_path)
        if len(df) == 0:
            print(f"DIGGER_WARNING: {locus}{sense_label} produced empty table (headers only): {csv_path}")
            return False
    except pd.errors.EmptyDataError:
        print(f"DIGGER_WARNING: {locus}{sense_label} produced completely empty file: {csv_path}")
        return False
    return True


def run_digger_per_locus(
    locus_fastas: dict[str, str],
    loci: dict[str, dict[str, str]],
    species: str,
    outdir: str,
    motif_dir: str | None = None,
) -> dict[str, str]:
    """Run digger on each per-locus FASTA file.

    IGK is handled specially: digger is run in both + and - sense
    (because IGK is inverted and duplicated) and the two output CSVs
    are merged into a single table.

    Returns a dict of locus -> path to the final digger output CSV.
    """
    digger_outdir = os.path.join(outdir, "digger_results")
    os.makedirs(digger_outdir, exist_ok=True)

    digger_outputs: dict[str, str] = {}

    for locus, fasta_path in locus_fastas.items():
        refs = loci[locus]
        locus_outdir = os.path.join(digger_outdir, locus)
        os.makedirs(locus_outdir, exist_ok=True)

        try:
            if locus == "IGK":
                # IGK is inverted and duplicated — run in both sense directions
                output_fwd = os.path.join(locus_outdir, f"{locus}_digger_output_fwd.csv")
                output_rev = os.path.join(locus_outdir, f"{locus}_digger_output_rev.csv")

                cmd_fwd = _build_digger_cmd(fasta_path, output_fwd, refs, species, locus, sense="+", motif_dir=motif_dir)
                print(f"Running digger for locus {locus} (sense +): {' '.join(cmd_fwd)}")
                subprocess.run(cmd_fwd, check=True, cwd=locus_outdir)
                _check_digger_output(output_fwd, locus, sense="+")

                cmd_rev = _build_digger_cmd(fasta_path, output_rev, refs, species, locus, sense="-", motif_dir=motif_dir)
                print(f"Running digger for locus {locus} (sense -): {' '.join(cmd_rev)}")
                subprocess.run(cmd_rev, check=True, cwd=locus_outdir)
                _check_digger_output(output_rev, locus, sense="-")

                # Merge the two IGK outputs
                merged_output = os.path.join(locus_outdir, f"{locus}_digger_output.csv")
                df_fwd = pd.read_csv(output_fwd)
                df_rev = pd.read_csv(output_rev)
                df_merged = pd.concat([df_fwd, df_rev], ignore_index=True)
                df_merged.to_csv(merged_output, index=False)
                print(f"Merged IGK forward/reverse digger outputs -> {merged_output}")
                digger_outputs[locus] = merged_output
            else:
                output_file = os.path.join(locus_outdir, f"{locus}_digger_output.csv")
                cmd = _build_digger_cmd(fasta_path, output_file, refs, species, locus, motif_dir=motif_dir)
                print(f"Running digger for locus {locus}: {' '.join(cmd)}")
                subprocess.run(cmd, check=True, cwd=locus_outdir)
                _check_digger_output(output_file, locus)
                digger_outputs[locus] = output_file
        except Exception as e:
            print(f"[WARNING] Processing failed for locus {locus}. Error: {e}")

    return digger_outputs


def map_reads_to_contigs(
    reads_fasta: str,
    contigs_fasta: str,
    outdir: str,
    minimap_option: str,
    threads: int,
) -> str:
    """Map CCS reads to the assembly contigs using minimap2.

    This mirrors the pers_to_ref mapping from get_read_support_VDJs.py:
    the contigs FASTA acts as the reference, CCS reads are mapped to it.

    Returns the path to the sorted, indexed BAM.
    """
    map_dir = os.path.join(outdir, "digger_read_support")
    os.makedirs(map_dir, exist_ok=True)

    sam_out = os.path.join(map_dir, "ccs_to_contigs.sam")
    bam_out = os.path.join(map_dir, "ccs_to_contigs.bam")
    sorted_bam = os.path.join(map_dir, "ccs_to_contigs.sorted.bam")

    # Index the contigs FASTA
    subprocess.run(["samtools", "faidx", contigs_fasta], check=True)

    # Map reads
    subprocess.run(
        [
            "minimap2", "-ax", minimap_option,
            "--secondary=yes", "-t", str(threads), "-L",
            contigs_fasta, reads_fasta,
            "-o", sam_out,
        ],
        check=True,
    )
    subprocess.run(["samtools", "view", "-Sbh", sam_out, "-o", bam_out], check=True)
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), bam_out, "-o", sorted_bam],
        check=True,
    )
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    # Clean up intermediates
    for f in [sam_out, bam_out]:
        try:
            os.remove(f)
        except OSError:
            pass

    return sorted_bam


def run_read_support(
    digger_outputs: dict[str, str],
    locus_fastas: dict[str, str],
    sorted_bam: str,
    outdir: str,
) -> None:
    """Run wasptk readsupport on each digger output CSV.

    wasptk readsupport -f <reference.fa> <allele_annotation.csv> <mapped.bam> <output.csv>
        --gene-col 'gene type' --start-col gene_start --end-col gene_end

    The reference FASTA is the per-locus contigs FASTA (i.e. the assembly
    contigs that digger annotated and that reads were mapped to).
    """
    rs_outdir = os.path.join(outdir, "digger_read_support")
    os.makedirs(rs_outdir, exist_ok=True)

    for locus, digger_csv in digger_outputs.items():
        if locus not in locus_fastas:
            print(f"[WARNING] No contigs FASTA for locus {locus}, skipping read support.")
            continue

        contigs_fasta = locus_fastas[locus]
        output_csv = os.path.join(rs_outdir, f"{locus}_digger_read_support.csv")

        cmd = [
            "wasptk", "readsupport",
            "-f", contigs_fasta,
            digger_csv,
            sorted_bam,
            output_csv,
            "--gene-col", "gene type",
            "--start-col", "start",
            "--end-col", "end",
        ]

        print(f"Running wasptk readsupport for {locus}: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Bin assembly contigs by locus, run digger per locus, and add read support."
    )
    parser.add_argument("outdir", help="Output directory (contains full_asm_for_digger.fasta)")
    parser.add_argument("-species", required=True, help="Species name matching reference file prefix (e.g. Homo_sapiens)")
    parser.add_argument("-allele_ref_dir", required=True, help="Directory containing species-prefixed allele reference FASTAs")
    parser.add_argument("-reads", required=True, help="Path to CCS reads FASTA for read support mapping")
    parser.add_argument("-minimap_option", default="map-hifi", help="Minimap2 preset (default: map-hifi)")
    parser.add_argument("-threads", type=int, default=4, help="Number of threads for minimap2/samtools")
    parser.add_argument("--locus_fasta", nargs="+", help="User-supplied contig fasta per locus. Format: LOCUS=FILE (e.g. IGH=path/to/igh.fasta)")
    parser.add_argument("--no-blast", action="store_true", help="Skip BLAST and only run digger on user-supplied locus fastas")
    parser.add_argument("-motif_dir", default=None, help="Path to motif directory for unsupported species")
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    species = args.species
    allele_ref_dir = os.path.abspath(args.allele_ref_dir)
    reads_fasta = os.path.abspath(args.reads)
    minimap_option = args.minimap_option
    threads = args.threads
    assembly_fasta = os.path.join(outdir, "full_asm_for_digger.fasta")

    user_fastas = {}
    if args.locus_fasta:
        for item in args.locus_fasta:
            if "=" not in item:
                sys.exit(f"ERROR: Invalid format for --locus_fasta: {item}. Expected LOCUS=FILE")
            loc, path = item.split("=", 1)
            if not os.path.isfile(path):
                sys.exit(f"ERROR: User-supplied fasta not found: {path}")
            user_fastas[loc] = os.path.abspath(path)

    print(f"Discovering loci from {allele_ref_dir} for species {species}...")
    loci = discover_loci(allele_ref_dir, species)
    print(f"Found loci: {', '.join(sorted(loci.keys()))}")

    loci_to_bin = {k: v for k, v in loci.items() if k not in user_fastas}
    locus_fastas = {}

    if args.no_blast:
        print("Skipping BLAST against assembly because --no-blast was specified.")
        loci_to_bin = {}  # Empty it so the binning step is skipped

    if loci_to_bin:
        if not os.path.isfile(assembly_fasta):
            sys.exit(f"ERROR: Assembly FASTA not found: {assembly_fasta}")
        print(f"Loci to bin from assembly: {', '.join(sorted(loci_to_bin.keys()))}")
        print("Building BLAST database from reference sequences...")
        db_path = build_blast_db(loci_to_bin, outdir)

        print("BLASTing contigs against reference database...")
        blast_results = blast_contigs(assembly_fasta, db_path, outdir)

        print("Binning contigs by locus...")
        locus_fastas = bin_contigs_by_locus(blast_results, assembly_fasta, loci_to_bin, outdir)
        print(f"Binned contigs into {len(locus_fastas)} loci: {', '.join(sorted(locus_fastas.keys()))}")

        # Warn about discovered loci that got zero contigs assigned
        missing_loci = set(loci_to_bin.keys()) - set(locus_fastas.keys())
        if missing_loci:
            for ml in sorted(missing_loci):
                print(f"BINNING_WARNING: No contigs were assigned to locus {ml}")
        if not locus_fastas:
            print("BLAST_WARNING: BLAST produced hits but no contigs could be assigned to any locus")

    # Add user supplied fastas
    for loc, path in user_fastas.items():
        if loc not in loci:
            print(f"[WARNING] User provided fasta for {loc} but no references found for this locus.")
        locus_fastas[loc] = path

    print("Running digger per locus...")
    digger_outputs = run_digger_per_locus(locus_fastas, loci, species, outdir, motif_dir=args.motif_dir)

    print("Mapping CCS reads to assembly contigs...")
    sorted_bam = map_reads_to_contigs(
        reads_fasta, assembly_fasta, outdir, minimap_option, threads
    )

    print("Running wasptk readsupport per locus...")
    run_read_support(digger_outputs, locus_fastas, sorted_bam, outdir)

    print("Done.")


if __name__ == "__main__":
    main()
