#!/usr/bin/env python3
"""
run_digger.py

Bins assembly contigs by IG/TR locus using BLAST against reference allele
sequences, then runs the digger tool on each per-locus FASTA.

Workflow:
  1. Scan allele_ref_dir for reference FASTAs matching the species prefix
     to discover which loci exist (e.g. IGH, IGK, IGL, TRA, TRB, ...).
  2. Build a combined BLAST database from all non-gapped reference FASTAs.
  3. BLAST the assembly contigs (full_asm_for_digger.fasta) against the DB.
  4. Assign each contig to its best-hit locus.
  5. Write per-locus FASTA files (e.g. IGH_contigs.fasta).
  6. Run digger on each per-locus FASTA, passing the appropriate
     species-specific V/D/J/C reference files.
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict


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


def build_blast_db(allele_ref_dir: str, species: str, outdir: str) -> str:
    """Concatenate all non-gapped reference FASTAs and build a BLAST DB.

    Returns the path to the BLAST database.
    """
    db_dir = os.path.join(outdir, "digger_blast_db")
    os.makedirs(db_dir, exist_ok=True)

    combined_fasta = os.path.join(db_dir, "combined_refs.fasta")
    pattern = os.path.join(allele_ref_dir, f"{species}_*.fasta")
    ref_files = sorted(
        f for f in glob.glob(pattern) if "_gapped" not in os.path.basename(f)
    )

    with open(combined_fasta, "w") as out_fh:
        for ref_file in ref_files:
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
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
            "-evalue", "1e-5",
            "-max_target_seqs", "1",
            "-out", blast_out,
        ],
        check=True,
    )
    return blast_out


def bin_contigs_by_locus(
    blast_results: str, assembly_fasta: str, loci: dict, outdir: str
) -> dict[str, str]:
    """Parse BLAST results, assign each contig to a locus, write per-locus FASTAs.

    Returns a dict of locus -> path to the per-locus contig FASTA.
    """
    # Map each contig to its best-hit locus based on the subject ID.
    # Subject IDs come from the reference FASTAs and contain the locus name,
    # e.g. a hit to "IGHV1-2*01" means the contig belongs to IGH.
    known_loci = set(loci.keys())

    contig_to_locus: dict[str, str] = {}
    with open(blast_results) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            qseqid, sseqid = fields[0], fields[1]
            # Try to extract the locus prefix from the subject sequence ID.
            # Reference sequence names typically start with the locus,
            # e.g. "IGHV1-2*01", "IGKJ5*01", "TRAV1*01"
            for locus in known_loci:
                if locus in sseqid:
                    contig_to_locus[qseqid] = locus
                    break

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

    for contig_id, locus in contig_to_locus.items():
        if contig_id in contigs:
            locus_contigs[locus].append(contigs[contig_id])

    for locus, seqs in locus_contigs.items():
        fasta_path = os.path.join(binned_dir, f"{locus}_contigs.fasta")
        with open(fasta_path, "w") as fh:
            fh.write("\n".join(seqs) + "\n")
        locus_fastas[locus] = fasta_path

    return locus_fastas


def run_digger_per_locus(
    locus_fastas: dict[str, str],
    loci: dict[str, dict[str, str]],
    species: str,
    outdir: str,
) -> None:
    """Run digger on each per-locus FASTA file."""
    digger_outdir = os.path.join(outdir, "digger_results")
    os.makedirs(digger_outdir, exist_ok=True)

    for locus, fasta_path in locus_fastas.items():
        refs = loci[locus]
        locus_outdir = os.path.join(digger_outdir, locus)
        os.makedirs(locus_outdir, exist_ok=True)

        output_file = os.path.join(locus_outdir, f"{locus}_digger_output.csv")

        cmd = [
            "digger",
            fasta_path,
            output_file,
            "-species", species,
        ]

        # Add reference file arguments if present
        if "V" in refs:
            cmd.extend(["-v_ref", refs["V"]])
        if "D" in refs:
            cmd.extend(["-d_ref", refs["D"]])
        if "J" in refs:
            cmd.extend(["-j_ref", refs["J"]])
        if "V_gapped" in refs:
            cmd.extend(["-v_ref_gapped", refs["V_gapped"]])

        # Determine the locus flag for digger (e.g. IGH, IGK, IGL)
        cmd.extend(["-locus", locus])

        print(f"Running digger for locus {locus}: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Bin assembly contigs by locus and run digger per locus."
    )
    parser.add_argument("outdir", help="Output directory (contains full_asm_for_digger.fasta)")
    parser.add_argument("-species", required=True, help="Species name matching reference file prefix (e.g. Homo_sapiens)")
    parser.add_argument("-allele_ref_dir", required=True, help="Directory containing species-prefixed allele reference FASTAs")
    args = parser.parse_args()

    outdir = args.outdir
    species = args.species
    allele_ref_dir = args.allele_ref_dir
    assembly_fasta = os.path.join(outdir, "full_asm_for_digger.fasta")

    if not os.path.isfile(assembly_fasta):
        sys.exit(f"ERROR: Assembly FASTA not found: {assembly_fasta}")

    print(f"Discovering loci from {allele_ref_dir} for species {species}...")
    loci = discover_loci(allele_ref_dir, species)
    print(f"Found loci: {', '.join(sorted(loci.keys()))}")

    print("Building BLAST database from reference sequences...")
    db_path = build_blast_db(allele_ref_dir, species, outdir)

    print("BLASTing contigs against reference database...")
    blast_results = blast_contigs(assembly_fasta, db_path, outdir)

    print("Binning contigs by locus...")
    locus_fastas = bin_contigs_by_locus(blast_results, assembly_fasta, loci, outdir)
    print(f"Binned contigs into {len(locus_fastas)} loci: {', '.join(sorted(locus_fastas.keys()))}")

    print("Running digger per locus...")
    run_digger_per_locus(locus_fastas, loci, species, outdir)

    print("Done.")


if __name__ == "__main__":
    main()
