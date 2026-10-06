#!/usr/bin/env python3
"""
cut_contigs.py

Trims assembly contigs that extend beyond IG/TR loci boundaries.

Given a mapped assembly BAM and a BED file defining IG/TR loci regions:
1. Expands each BED interval by a buffer distance (default 20kb)
2. Merges overlapping expanded intervals
3. For each aligned contig, trims portions extending beyond the merged intervals
4. Writes trimmed contigs to a new FASTA

Usage:
    python cut_contigs.py <bam> <bed> <output_fasta> [--buffer 20000]
"""

import argparse
import sys
from collections import defaultdict

import pysam


def read_bed(bed_path):
    """Read a BED file and return a list of (chrom, start, end, name) tuples."""
    intervals = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else ""
            intervals.append((chrom, start, end, name))
    return intervals


def expand_and_merge(intervals, buffer_distance):
    """Expand intervals by buffer, then merge overlapping ones per chromosome.

    Returns a dict of chrom -> list of (start, end) merged intervals.
    """
    by_chrom = defaultdict(list)
    for chrom, start, end, _name in intervals:
        expanded_start = max(0, start - buffer_distance)
        expanded_end = end + buffer_distance
        by_chrom[chrom].append((expanded_start, expanded_end))

    merged = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        merged_list = [list(ivs[0])]
        for start, end in ivs[1:]:
            if start <= merged_list[-1][1]:
                merged_list[-1][1] = max(merged_list[-1][1], end)
            else:
                merged_list.append([start, end])
        merged[chrom] = [(s, e) for s, e in merged_list]

    return merged


def find_cut_boundaries(ref_start, ref_end, chrom_intervals):
    """Determine reference-coordinate cut boundaries for a contig.

    Returns (left_cut_ref, right_cut_ref).  None means no cut on that side.
    """
    if not chrom_intervals:
        return None, None

    overlapping = [
        (s, e) for s, e in chrom_intervals if s < ref_end and e > ref_start
    ]

    if not overlapping:
        return None, None

    allowed_start = min(s for s, _ in overlapping)
    allowed_end = max(e for _, e in overlapping)

    left_cut = allowed_start if ref_start < allowed_start else None
    right_cut = allowed_end if ref_end > allowed_end else None

    return left_cut, right_cut


def ref_to_query_pos(aligned_pairs, ref_pos, side="left"):
    """Translate a reference position to a query position via aligned pairs.

    For 'left' cuts: returns the first query pos at or after ref_pos.
    For 'right' cuts: returns the last query pos at or before ref_pos.
    """
    if not aligned_pairs:
        return None

    if side == "left":
        for qpos, rpos in aligned_pairs:
            if rpos >= ref_pos:
                return qpos
        return None
    else:
        result = None
        for qpos, rpos in aligned_pairs:
            if rpos <= ref_pos:
                result = qpos
            else:
                break
        return result


def main():
    parser = argparse.ArgumentParser(
        description="Trim assembly contigs extending beyond IG/TR loci."
    )
    parser.add_argument("bam", help="Mapped assembly BAM")
    parser.add_argument("bed", help="BED file defining IG/TR loci regions")
    parser.add_argument("output_fasta", help="Output FASTA with trimmed contigs")
    parser.add_argument(
        "--buffer", type=int, default=20000,
        help="Buffer distance to extend BED regions in bp (default: 20000)"
    )
    args = parser.parse_args()

    intervals = read_bed(args.bed)
    merged = expand_and_merge(intervals, args.buffer)

    print(f"Loaded {len(intervals)} BED intervals, expanded by {args.buffer}bp and merged:")
    for chrom in sorted(merged):
        for start, end in merged[chrom]:
            print(f"  {chrom}:{start}-{end}")

    bam = pysam.AlignmentFile(args.bam, "rb")

    written = 0
    skipped = 0
    cut_count = 0
    unchanged = 0

    with open(args.output_fasta, "w") as out_fh:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            contig_name = read.query_name
            seq = read.query_sequence
            if seq is None:
                continue

            ref_name = read.reference_name
            ref_start = read.reference_start
            ref_end = read.reference_end

            chrom_ivs = merged.get(ref_name, [])
            left_cut_ref, right_cut_ref = find_cut_boundaries(
                ref_start, ref_end, chrom_ivs
            )

            if left_cut_ref is None and right_cut_ref is None:
                out_fh.write(f">{contig_name}\n{seq}\n")
                written += 1
                unchanged += 1
                continue

            aligned_pairs = [
                (qpos, rpos)
                for qpos, rpos in read.get_aligned_pairs(matches_only=True)
            ]

            if not aligned_pairs:
                out_fh.write(f">{contig_name}\n{seq}\n")
                written += 1
                unchanged += 1
                continue

            q_start = 0
            q_end = len(seq)

            if left_cut_ref is not None:
                qpos = ref_to_query_pos(aligned_pairs, left_cut_ref, side="left")
                if qpos is not None:
                    q_start = qpos

            if right_cut_ref is not None:
                qpos = ref_to_query_pos(aligned_pairs, right_cut_ref, side="right")
                if qpos is not None:
                    q_end = qpos + 1

            trimmed_seq = seq[q_start:q_end]

            if len(trimmed_seq) == 0:
                print(f"CUT_WARNING: {contig_name} was trimmed to 0bp, skipping")
                skipped += 1
                continue

            original_len = len(seq)
            new_len = len(trimmed_seq)
            if new_len < original_len:
                trimmed_bp = original_len - new_len
                print(
                    f"CUT_INFO: {contig_name} trimmed from {original_len} to "
                    f"{new_len} ({trimmed_bp}bp removed) "
                    f"(ref {ref_name}:{ref_start}-{ref_end})"
                )
                cut_count += 1

            out_fh.write(f">{contig_name}\n{trimmed_seq}\n")
            written += 1

    bam.close()

    print(
        f"\nDone. Wrote {written} contigs: "
        f"{cut_count} trimmed, {unchanged} unchanged, {skipped} skipped."
    )


if __name__ == "__main__":
    main()
