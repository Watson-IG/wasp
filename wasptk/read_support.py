"""Compute read support metrics from a mapped BAM."""

from __future__ import annotations

from typing import List, Dict, Tuple, Set

import pandas as pd
import pysam

from .match_subsequences import extract_sequence, count_matching_reads as count_match_sub
from .ighc_match import extract_exon_sequences, count_matching_reads as count_ighc

_CONSTANT_GENES = ["IGKC", "IGLC", "TRAC", "TRBC", "TRDC", "TRGC", "IGHC"]


def _mismatch_positions(read: pysam.AlignedSegment) -> Set[int]:
    """Return reference positions with mismatches using the MD tag."""
    try:
        md = read.get_tag("MD")
    except KeyError:
        return set()
    pos = read.reference_start
    mismatches: Set[int] = set()
    i = 0
    while i < len(md):
        if md[i].isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            pos += int(md[i:j])
            i = j
        elif md[i] == "^":
            j = i + 1
            while j < len(md) and md[j].isalpha():
                j += 1
            pos += j - i - 1
            i = j
        else:
            # one or more mismatch letters
            while i < len(md) and md[i].isalpha() and md[i] != "^":
                mismatches.add(pos)
                pos += 1
                i += 1
    return mismatches


def _pileup_region(
    bam: pysam.AlignmentFile, contig: str, start: int, end: int
) -> Tuple[List[int], List[int], List[int]]:
    """Return coverage, mismatch counts and match counts for a region."""
    length = end - start + 1
    coverage = [0] * length
    mismatches = [0] * length
    matches = [0] * length
    for read in bam.fetch(contig, start - 1, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        mismatch_pos = _mismatch_positions(read)
        for pos in read.get_reference_positions():
            if start - 1 <= pos <= end - 1:
                idx = pos - (start - 1)
                coverage[idx] += 1
                if pos in mismatch_pos:
                    mismatches[idx] += 1
                else:
                    matches[idx] += 1
    return coverage, mismatches, matches


def _calc_simple(
    coverage: List[int], mism: List[int], match: List[int]
) -> Tuple[int, int, float, int, int, str, str, float, int]:
    total_positions = len(coverage)
    total_reads = sum(coverage)
    mismatched_positions = 0
    matched_positions = 0
    positions_with_10x = 0
    mismatch_list: List[str] = []
    match_list: List[str] = []
    for c, m, ma in zip(coverage, mism, match):
        mismatch_list.append(str(m))
        match_list.append(str(ma))
        if c >= 10:
            positions_with_10x += 1
        mismatch_rate = m / c if c > 0 else 0
        match_rate = ma / c if c > 0 else 0
        if mismatch_rate > 0.2:
            mismatched_positions += 1
        if match_rate > 0.8:
            matched_positions += 1
    avg_reads_per_position = total_reads / total_positions if total_positions else 0
    percent_accuracy = (
        (matched_positions / total_positions) * 100 if total_positions else 0
    )
    return (
        total_positions,
        total_reads,
        avg_reads_per_position,
        mismatched_positions,
        matched_positions,
        ":".join(mismatch_list),
        ":".join(match_list),
        percent_accuracy,
        positions_with_10x,
    )


def _calc_ighc(
    coverage: List[int], mism: List[int], match: List[int]
) -> Tuple[int, int, float, int, int, str, str, int, int, int, int, int, float]:
    total_positions = len(coverage)
    total_reads = sum(coverage)
    mismatched_positions = 0
    matched_positions = 0
    mismatch_list: List[str] = []
    match_list: List[str] = []
    mismatched_positions_cov_lt10 = 0
    mismatched_positions_cov_ge10 = 0
    matched_positions_cov_lt10 = 0
    matched_positions_cov_ge10 = 0
    pos_10x = 0
    for c, m, ma in zip(coverage, mism, match):
        mismatch_list.append(str(m))
        match_list.append(str(ma))
        if c >= 10:
            pos_10x += 1
        mismatch_rate = m / c if c > 0 else 0
        match_rate = ma / c if c > 0 else 0
        if mismatch_rate > 0.2:
            mismatched_positions += 1
            if c < 10:
                mismatched_positions_cov_lt10 += 1
            else:
                mismatched_positions_cov_ge10 += 1
        if match_rate >= 0.8:
            matched_positions += 1
            if c < 10:
                matched_positions_cov_lt10 += 1
            else:
                matched_positions_cov_ge10 += 1
    percent_accuracy = (
        (matched_positions / total_positions) * 100 if total_positions else 0
    )
    avg_cov = total_reads / total_positions if total_positions else 0
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
        percent_accuracy,
    )


def compute_read_support(
    allele_table: str,
    bam_path: str,
    output: str,
    contig_col: str = "contig",
    start_col: str = "start",
    end_col: str = "end",
    gene_col: str = "gene",
) -> None:
    df = pd.read_csv(allele_table)
    bam = pysam.AlignmentFile(bam_path, "rb")

    metrics: List[Dict[str, int | float | str]] = []

    for _, row in df.iterrows():
        contig = str(row[contig_col])
        start = int(row[start_col])
        end = int(row[end_col])
        gene = str(row.get(gene_col, ""))

        if any(gene.startswith(c) for c in _CONSTANT_GENES):
            seqs, coords = extract_exon_sequences(row)
            cov_all: List[int] = []
            mism_all: List[int] = []
            match_all: List[int] = []
            if coords:
                for s, e in coords:
                    cov, mism, match = _pileup_region(bam, contig, s, e)
                    cov_all.extend(cov)
                    mism_all.extend(mism)
                    match_all.extend(match)
            else:
                cov_all, mism_all, match_all = _pileup_region(bam, contig, start, end)
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
                pct_acc,
            ) = _calc_ighc(cov_all, mism_all, match_all)
            full_span, all_match, match_counts, span_counts = count_ighc(
                bam, contig, start, end, seqs, coords
            )
            row_metrics: Dict[str, int | float | str] = {
                "Total_Positions": tpos,
                "Total_Reads_by_Positions": total_reads,
                "Average_Coverage": avg_cov,
                "Mismatched_Positions": mismatched_positions,
                "Matched_Positions": matched_positions,
                "Position_Mismatches": mismatch_str,
                "Position_Matches": match_str,
                "Mismatched_Positions_Coverage_Less_Than_10": mismatched_positions_cov_lt10,
                "Mismatched_Positions_Coverage_10_Or_Greater": mismatched_positions_cov_ge10,
                "Matched_Positions_Coverage_Less_Than_10": matched_positions_cov_lt10,
                "Matched_Positions_Coverage_10_Or_Greater": matched_positions_cov_ge10,
                "Positions_With_At_Least_10x_Coverage": pos_10x,
                "Percent_Accuracy": pct_acc,
                "Fully_Spanning_Reads": full_span,
                "Fully_Spanning_Reads_100%_Match": all_match,
            }
            for i, count in enumerate(match_counts, start=1):
                row_metrics[f"Allele_reads_100_Match_e{i}"] = count
            for i, count in enumerate(span_counts, start=1):
                row_metrics[f"Allele_reads_fully_spanning_e{i}"] = count
        else:
            cov, mism, match = _pileup_region(bam, contig, start, end)
            (
                tpos,
                total_reads,
                avg_cov,
                mismatched_positions,
                matched_positions,
                mismatch_str,
                match_str,
                pct_acc,
                pos_10x,
            ) = _calc_simple(cov, mism, match)
            seq = extract_sequence(row, gene)
            full_span, perfect = count_match_sub(bam, contig, start, end, seq)
            row_metrics = {
                "Total_Positions": tpos,
                "Total_Reads_by_Positions": total_reads,
                "Average_Coverage": avg_cov,
                "Mismatched_Positions": mismatched_positions,
                "Matched_Positions": matched_positions,
                "Position_Mismatches": mismatch_str,
                "Position_Matches": match_str,
                "Percent_Accuracy": pct_acc,
                "Positions_With_At_Least_10x_Coverage": pos_10x,
                "Fully_Spanning_Reads": full_span,
                "Fully_Spanning_Reads_100%_Match": perfect,
            }

        metrics.append(row_metrics)

    all_cols: List[str] = []
    for m in metrics:
        for k in m.keys():
            if k not in all_cols:
                all_cols.append(k)

    for col in all_cols:
        df[col] = [m.get(col, 0) for m in metrics]

    df.to_csv(output, index=False)
