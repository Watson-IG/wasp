#!/usr/bin/env python3

import csv
import sys
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])

def get_sequences_and_regions_from_csv(import_csv, gene_key, contig):
    sequences = []
    regions = []
    exon_start_end_list = []
    with open(import_csv, mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if row['gene'] == gene_key and row['contig'] == contig:
                start = int(float(row['allele_sequence_start']))
                end = int(float(row['allele_sequence_end']))
                reverse_comp = 'sense' in row and row['sense'] == '-'
                # Populate sequences from C-EXON columns
                for i in range(1, 10):
                    sequence_column = f'C-EXON_{i}'
                    seq_start = f'C-EXON_{i}_start'
                    seq_end = f'C-EXON_{i}_end'
                    sequence = row.get(sequence_column, '')
                    if row.get(seq_start, '') != '' and row.get(seq_end, '') != '':
                        sequence_start_end = (int(float(row.get(seq_start, ''))), int(float(row.get(seq_end, '')))) # exon specific start-end tuples
                    if reverse_comp and sequence:
                        sequence = reverse_complement(sequence)
                    if sequence:  # Only add non-empty sequences
                        sequences.append(sequence)
                        exon_start_end_list.append(sequence_start_end)
                regions.append((contig, start, end))  # Make sure to append a tuple of (contig, start, end)
                break
    return sequences, regions, exon_start_end_list

def count_matching_reads(bamfile, sequences, regions, exon_start_end_list):
    full_span_counts = []
    full_span_all_match_count_list = []
    perfect_match_counts = []
    perfect_spans_counts = []
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for region in regions:
        start, end = region[1], region[2]  # Assuming region tuple structure is (contig, start, end)
        contig = str(region[0])  # Ensure this is a string as expected
        full_span_count = 0
        perfect_matches = [0] * len(sequences)
        perfect_spans = [0] * len(sequences)
        all_matches = [0] * len(sequences)
        full_span_all_match_count = 0
        for read in samfile.fetch(str(contig), start, end):  # Pass contig correctly
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_start <= start and read.reference_end >= end:
                full_span_count += 1
                read_seq = read.query_sequence
                for idx, sequence in enumerate(sequences):
                    if sequence and (sequence in read_seq or reverse_complement(sequence) in read_seq):
                        all_matches[idx] = 1
                if 0 not in all_matches:
                    full_span_all_match_count += 1
            read_seq = read.query_sequence
            for idx, start_end in enumerate(exon_start_end_list):
                ex_start = start_end[0]
                ex_end = start_end[1]
                if read.reference_start <= ex_start and read.reference_end >= ex_end:
                    perfect_spans[idx] += 1
            for idx, sequence in enumerate(sequences):
                if sequence and (sequence in read_seq or reverse_complement(sequence) in read_seq): 
                    perfect_matches[idx] += 1
        full_span_counts.append(full_span_count)
        full_span_all_match_count_list.append(full_span_all_match_count)
        perfect_match_counts.append(perfect_matches)
        perfect_spans_counts.append(perfect_spans)
    return full_span_counts, full_span_all_match_count_list, perfect_match_counts, perfect_spans_counts

if __name__ == '__main__':
    bamfile = sys.argv[1]
    contig = sys.argv[2]
    gene_key = sys.argv[3]
    import_csv = sys.argv[4]

    sequences, regions, exon_start_end_list = get_sequences_and_regions_from_csv(import_csv, gene_key, contig)
    if sequences and regions:
        full_span_counts, full_span_all_match_count_list, perfect_match_counts, perfect_spans_counts = count_matching_reads(bamfile, sequences, regions, exon_start_end_list) 
        for i, (i_counts, j_counts) in enumerate(zip(perfect_match_counts, perfect_spans_counts)):
            i_counts_extended = (i_counts + [""] * 9)[:9] # fill 9 exon pos
            j_counts_extended = (j_counts + [""] * 9)[:9]
            
            print(f"{full_span_counts[i]},{full_span_all_match_count_list[i]},{','.join(map(str, i_counts_extended))},{','.join(map(str, j_counts_extended))}") # first we will get 100% exon matches then #no of fully spanning exon reads

        #for i, counts in enumerate(perfect_match_counts):
        #    print(f"{full_span_counts[i]},{full_span_all_match_count[i]},{','.join(map(str, counts))}")
    else:
        print(f"No sequences or regions found for gene {gene_key} on contig {contig}", file=sys.stderr)
