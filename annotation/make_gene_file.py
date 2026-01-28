# Create a file of identified genes from a set of assemblies
import argparse
import csv
from collections import namedtuple

from read_bed import get_gene_type, read_beds
import pysam
import receptor_utils.simple_bio_seq as simple
import glob
import os

from read_indexed_fasta import read_sequence_from_fasta

from cigar import Cigar

output_headers = [
    'project',
    'subject',
    'sample_name',
    'haplotype',
    'gene',
    'allele',
]

v_coord_map = {
    'allele_sequence': 'GENE',
    'V-HEPTAMER': 'HEPTAMER',
    'V-SPACER': 'SPACER',
    'V-NONAMER': 'NONAMER',
    'V-EXON2': 'EXON_2',
    'V-INTRON': 'INTRON',
    'L-PART1': 'EXON_1',
    'V-REGION': 'V-REGION',
    'L-PART2': 'L-PART2',
    'V-UTR': 'UTR',
}

utr_loci = ['IGH']

d_coord_map = {
    'allele_sequence': 'GENE',
    'D-3_HEPTAMER': '3_HEPTAMER',
    'D-3_SPACER': '3_SPACER',
    'D-3_NONAMER': '3_NONAMER',
    'D-REGION': 'EXON_1',
    'D-5_HEPTAMER': '5_HEPTAMER',
    'D-5_SPACER': '5_SPACER',
    'D-5_NONAMER': '5_NONAMER',
}

j_coord_map = {
    'allele_sequence': 'GENE',
    'J-HEPTAMER': 'HEPTAMER',
    'J-SPACER': 'SPACER',
    'J-NONAMER': 'NONAMER',
    'J-REGION': 'EXON_1',
}

c_coord_map = {
}


def missing_fields(row, fields):
    for field in fields:
        if field not in row or not row[field] or len(row[field]) < 2:
            return True
    return False


def verbose_print(verbose_str, verbose):
    if verbose:
        print(verbose_str)


Contig = namedtuple('Contig', 'header, annotations, haplotype, coords')


class CigarState:
    def __init__(self):
        self.count = 0
        self.state = 'M'
        self.next_pos = 0
        self.string = ''
        self.processed_seq_len = 0
        self.this_seg = ''
        self.processed_segs = []

    def process_cigar(self, op, c):
        if self.state != op:
            if self.count > 0:
                self.string += f"{self.count}{self.state}"
                self.processed_seq_len += self.count
                self.processed_segs.append((len(self.this_seg), self.this_seg))
            self.count = 1
            self.this_seg = c if op != 'D' else ''
            self.state = op
        else:
            if op != 'D':
                self.this_seg += c
            self.count += 1
        if op != 'I':
            self.next_pos += 1

    def finalise_cigar(self):
        if self.count > 0:
            self.string += f"{self.count}{self.state}"
            self.processed_seq_len += self.count
            self.processed_segs.append((len(self.this_seg), self.this_seg))
        return self.string


HAPLOTYPE_WARNING = False


# Fetch all contigs from the BAM file that cover the specified range
def fetch_contigs(samfile, chrom, start, end, annotation_ranges, verbose, list_all_contigs):
    contigs = []
    dupes = {}

    # if verbose:
    #    breakpoint()

    for read in samfile.fetch(chrom, start, end):
        if read.is_secondary:
            verbose_print(f"{read.qname}: secondary alignment", verbose)
            continue
        if read.is_supplementary:
            verbose_print(f"{read.qname}: supplementary read", verbose)
            continue
        if read.is_unmapped:
            verbose_print(f"{read.qname}: unmapped read", verbose)
            continue
        if read.reference_start > start:
            verbose_print(f"{read.qname} start: {read.reference_start} end {read.reference_end} read.ref_start > required start", verbose)
            continue
        if read.reference_end < end:
            verbose_print(f"{read.qname} start: {read.reference_start} end {read.reference_end} read.ref_end < required end", verbose)
            continue

        q_start = -1
        q_end = None
        leading_dels = 0

        picked_pairs = []
        for query, ref in read.get_aligned_pairs():
            if ref and ref == start:    # q_start is where we want to take our query sequence from...
                q_start = query         # but if it is None at this point, the first base is a deletion
            if q_start is None and query is not None and ref is not None and ref > start:   # in which case, we set q_start to the first base that isn't deleted...
                q_start = query
            if ref is not None and q_start is None:                     # and in the meantime count the leading deletions
                leading_dels += 1
            if query is not None and (ref is None or ref <= end):       # make sure q_end never runs past the end, but always runs right up to it
                q_end = query

            if q_start and q_start >= 0:
                picked_pairs.append((query, ref))

            if ref and ref >= end:
                break

        name = read.query_name

        if q_start is None or q_end is None:
            continue    # we are not going to achieve anything if there is no read coverage at all over the required range
        else:
            seq = read.query_sequence[q_start:q_end+1]

        haplotype = None
        for el in name.split('_'):
            if 'h=' in el:
                haplotype = el.replace('h=', '')        # this is the IGenotyper convention

        if haplotype is None:                           # hap1, hap2 is the Wasp convention
            if 'hap_1' in name or 'hap1' in name:
                haplotype = '1'
            elif 'hap_2' in name or 'hap2' in name:
                haplotype = '2'

        if haplotype is None or haplotype not in ['0', '1', '2']:
            print(f"{read.qname}: unrecognisable haplotype", verbose)
            haplotype = '0'

        # assign to annotation fields, handling indels
        annotations = {annot: '' for _, _, annot in annotation_ranges}
        contig_coords = {annot: list() for _, _, annot in annotation_ranges}

        cigar_states = {annot: CigarState() for _, _, annot in annotation_ranges}
        ap = read.get_aligned_pairs()

        if ap is None or (ap[0][1] is not None and ap[0][1] >= end):
            verbose_print(f"{read.qname}: no useful alignment with ref", verbose)
            continue    # no useful alignment

        def add_to_annots(pos, c, op, q_pos):
            for r in annotation_ranges:
                if pos >= r[0] and pos <= r[1]:
                    annotations[r[2]] += c
                    cigar = cigar_states[r[2]]
                    cigar.process_cigar(op, c)
                    if len(contig_coords[r[2]]) == 0:
                        contig_coords[r[2]].append(q_pos)
                        contig_coords[r[2]].append(q_pos)
                    elif len(contig_coords[r[2]]) == 2:
                        contig_coords[r[2]][1] = q_pos               

        def finalise_cigars():
            for r in annotation_ranges:
                cigar = cigar_states[r[2]]
                annotations[r[2] + '_CIGAR'] = cigar.finalise_cigar()
                cig = Cigar(annotations[r[2] + '_CIGAR'])
                if cig.__len__() != len(annotations[r[2]]):
                    print(f'Error in cigar string length for annotation {r[2]}')

        # bump up each coord by 1 to make them 1-based
        def finalise_coords(contig_coords):
            for k, v in contig_coords.items():
                if len(v) == 2 and v[0] is not None and v[1] is not None:
                    contig_coords[k] = (v[0] + 1, v[1] + 1)

        q_current = q_start
        seq_current = 0
        annot_pos = 1   # 1-based index of gene sequence relative to ref
        start_encountered = False
        wanted_length = end - start

        for _ in range(leading_dels):
            add_to_annots(annot_pos, '', 'D', q_current)
            annot_pos += 1

        for i in range(len(ap)):
            if start_encountered and (ap[i][0] is None or seq_current >= len(seq)):  # deletion relative to franken (including a trailing deletion)
                add_to_annots(annot_pos, '', 'D', q_current)
                annot_pos += 1
            elif ap[i][1] is None and ap[i][0] is not None and start_encountered:  # insertion relative to franken
                add_to_annots(annot_pos, seq[seq_current], 'I', q_current)
                seq_current += 1
                q_current += 1
            elif ap[i][0] and ap[i][0] == q_current:    # aligned relative to franken
                start_encountered = True
                add_to_annots(annot_pos, seq[seq_current], 'M', q_current)
                seq_current += 1
                q_current += 1
                annot_pos += 1

            if annot_pos > wanted_length:
                break

        finalise_cigars()
        verbose_print(f"{read.qname}: h={haplotype} {seq}", verbose)
        new_contig = Contig(name, annotations, haplotype, contig_coords)
        finalise_coords(new_contig.coords)

        contigs.append(new_contig)

    # remove hap 0 contigs if we have some that are hap 1 or hap 2

    non_zero_seen = False

    for contig in contigs:
        if contig.haplotype != '0':
            non_zero_seen = True
            break

    if non_zero_seen:
        for i in range(len(contigs)-1, -1, -1):
            if contigs[i].haplotype == '0':
                del contigs[i]

    # if there is at least one conting in a haplotype that has no Ds in the CIGAR, remove any contigs in that haplotype that do have Ds

    for hap in ['0', '1', '2']:
        clean_seen = False
        for contig in contigs:
            if contig.haplotype == hap and 'D' not in contig.annotations['allele_sequence_CIGAR']:
                clean_seen = True
                break
        if clean_seen:
            for i in range(len(contigs)-1, -1, -1):
                if contigs[i].haplotype == hap and 'D' in contigs[i].annotations['allele_sequence_CIGAR']:
                    del contigs[i]

    # dedupe unlesss we've been asked for everything

    dupes = {}
    if not list_all_contigs:
        final_contigs = list()
        for hap in ['0', '1', '2']:
            found_contig = None
            for contig in contigs:
                if contig.haplotype == hap:
                    if not found_contig:
                        found_contig = contig
                        dupes[found_contig.header] = list()
                    else:
                        dupes[found_contig.header].append(contig.header)
            
            if found_contig:
                final_contigs.append(found_contig)
        contigs = final_contigs

    if verbose:
        print(f"Contigs added to analysis {('excluding duplicates' if not list_all_contigs else 'including duplicates')}:")
        for contig in contigs:
            if 'v-exon2' in contig.annotations:
                print(f"{contig.header} {contig.haplotype} exon_2: {contig.annotations['v-exon2']}")

    return contigs, dupes


def process_rows(required_gene_type, refname, coord_map, samfile, project, subject, sample_name, beds, debug_gene, sense, forced_haplotype, list_all_contigs):
    rows = []

    for gene in beds[refname].keys():
        gene_type = get_gene_type(gene)

        if gene_type == 'C' and '_' in gene and gene[-2:] != '_1':
            continue        # for multi-exon genes, process all exons when we hit the first one

        if gene_type == required_gene_type:
            annotation_ranges = []

            rec_sense = sense
            if sense == '+-':
                if gene_type == 'V':
                    rec_sense = '+' if beds[refname][gene]['NONAMER']['start'] > beds[refname][gene]['HEPTAMER']['start'] else '-'
                elif gene_type == 'J':
                    rec_sense = '-' if beds[refname][gene]['NONAMER']['start'] > beds[refname][gene]['HEPTAMER']['start'] else '+' 
                elif gene_type == 'D':
                    for ex_gene in beds[refname]:       # find the first J gene and use its sense
                        if get_gene_type(ex_gene) == 'J':
                            rec_sense = '-' if beds[refname][ex_gene]['NONAMER']['start'] > beds[refname][ex_gene]['HEPTAMER']['start'] else '+'
                            break
                elif gene_type == 'C':
                    rec_sense = None    # decide below, based on order of exons

            if gene_type != 'C':
                seq_start = beds[refname][gene]['GENE']['start']
                seq_end = beds[refname][gene]['GENE']['end']
                gene_label = gene

                for annot_key, b_name in coord_map.items():
                    annotation_ranges.append((beds[refname][gene][b_name]['start'] - seq_start + 1, beds[refname][gene][b_name]['end'] - seq_start, annot_key))         # 1-based ranges
            else:
                seq_start = None
                seq_end = None
                rec_sense = None

                def get_gene_label(gene):
                    if '_' in gene:
                        return '_'.join(gene.split('_')[:-1])
                    else:
                        return gene

                gene_label = get_gene_label(gene)
                for ex_gene in beds[refname]:
                    if get_gene_type(ex_gene) == 'C' and gene_label == get_gene_label(ex_gene):
                        if seq_start is None or beds[refname][ex_gene]['GENE']['start'] < seq_start:
                            seq_start = beds[refname][ex_gene]['GENE']['start']

                        if seq_end is None or beds[refname][ex_gene]['GENE']['end'] > seq_end:
                            seq_end = beds[refname][ex_gene]['GENE']['end']

                        if rec_sense is None and '_' in ex_gene and ex_gene.split("_")[-1] != '1':
                            rec_sense = '+' if beds[refname][ex_gene]['GENE']['start'] > beds[refname][gene]['GENE']['start'] else '-'

                # For single exon C-genes, should there be any, find the first J gene and use its sense
                if rec_sense is None:
                    for ex_gene in beds[refname]:
                        if get_gene_type(ex_gene) == 'J':
                            rec_sense = '-' if beds[refname][ex_gene]['NONAMER']['start'] > beds[refname][ex_gene]['HEPTAMER']['start'] else '+'
                            break

                annotation_ranges.append((1, seq_end - seq_start, 'allele_sequence'))

                for ex_gene in beds[refname]:
                    if gene_label == '_'.join(ex_gene.split('_')[:-1]):
                        annotation_ranges.append((beds[refname][ex_gene]['GENE']['start'] - seq_start + 1, beds[refname][ex_gene]['GENE']['end'] - seq_start, f'C-EXON_{ex_gene.split("_")[-1]}'))

            if gene == debug_gene:
                print(f"Processing {gene} for {sample_name}. Required range {seq_start} - {seq_end}")

            contigs, dupes = fetch_contigs(samfile, refname, seq_start, seq_end, annotation_ranges, gene == debug_gene, list_all_contigs)

            # TODO - consider trying to fetch a contig for the coding region, if we can't get the whole gene
            
            for contig in contigs:
                contig_names = contig.header
                if contig.header in dupes:
                    contig_names += ',' + ','.join(dupes[contig.header])
                row = {
                    'project': project,
                    'subject': subject,
                    'sample_name': sample_name,
                    'contig': contig_names,
                    'haplotype': forced_haplotype if forced_haplotype else contig.haplotype,
                    'gene': gene_label,
                    'allele': '',
                    'sense': rec_sense,
                }
                for k, v in contig.annotations.items():
                    row[k] = v

                for k, v in contig.coords.items():
                    if len(v) == 2 and v[0] is not None and v[1] is not None:
                        row[k + '_start'] = v[0]
                        row[k + '_end'] = v[1]

                rows.append(row)

            # get contigs just covering the V-exon for comparison
            '''if label == 'HV':
                start =  beds['igh'][gene]['exon_2']['start']
                end = beds['igh'][gene]['exon_1']['end']
                exons = fetch_contigs(samfile, 'igh', seq_start, seq_end, [(1, 1 + end - start, 'allele_sequence')], list_all_contigs, gene == debug_gene)
                if len(contigs) != len(exons):
                    print(f"{gene}: full length {len(contigs)} exon only {len(exons)}") 
                    '''
    return rows


def process_sample(locus, refname, project, subject, sample_file, sample_name, beds, debug_gene, sense, forced_haplotype, list_all_contigs):
    samfile = pysam.AlignmentFile(sample_file)
    rows = []
    rows.extend(process_rows('V', refname, v_coord_map, samfile, project, subject, sample_name, beds, debug_gene, sense, forced_haplotype, list_all_contigs))

    if locus in ['IGH', 'TRB', 'TRD']:
        rows.extend(process_rows('D', refname, d_coord_map, samfile, project, subject, sample_name, beds, debug_gene, sense, forced_haplotype, list_all_contigs))

    rows.extend(process_rows('J', refname, j_coord_map, samfile, project, subject, sample_name, beds, debug_gene, sense, forced_haplotype, list_all_contigs))
    rows.extend(process_rows('C', refname, c_coord_map, samfile, project, subject, sample_name, beds, debug_gene, sense, forced_haplotype, list_all_contigs))

    for row in rows:
        for el in output_headers:
            if el not in row:
                row[el] = ''

    return rows

# sense defines the sense of the alleles within the reference sequence. + means 5' to 3' order.
# '+-' is used for loci that contain both + and - sense alleles, eg IGK. In this case, the sense of V and J genes will
# be determined from the position of the rss relative to the coding region. For other genes, the + or - sense must be specified
# in the 5th (last) column of the bed file containing REGION co-ordinates for the gene 
def main():
    parser = argparse.ArgumentParser(description='Extract VDJ gene fields from BAM file(s)')
    parser.add_argument('locus', help='locus (e.g. IGH, IGL)')
    parser.add_argument('refname', help='reference assembly name in SAM files (e.g. chr22)')
    parser.add_argument('sense', help='sense of the genes in the reference sequence: - is 5 to 3, + is 3 to 5, +- is both (use this for K only)')
    parser.add_argument('bed_dir', help='pathname to a directory holding one or more BED files containing the field co-ordinates')
    parser.add_argument('ref_seq', help='pathname to a FASTA file holding the reference sequence (in the same orientation as the BED files use it)')
    parser.add_argument('outfile', help='output file that will contain all the VDJ gene field data (csv)')
    parser.add_argument('--debug_gene', help='print debug output on the processing of this gene')
    parser.add_argument('--bam', '-b', help='pathname to a single BAM file')
    parser.add_argument('--project', '-p', default='Unknown', help='project name for the single BAM file')
    parser.add_argument('--subject', '-s', default='Unknown', help='subject name for the single BAM file')
    parser.add_argument('--sample', '-n', default='Unknown', help='sample name for the single BAM file')
    parser.add_argument('--bam_dirs', '-d', help='pathname to a directory BAM files (structure is project/subject/sample)')
    parser.add_argument('--list_all_contigs', '-l', help='list all contigs covering each gene on separate rows rather than one row per haplotype', action='store_true')
    args = parser.parse_args()

    locus = args.locus
    refname = args.refname
    sense = args.sense

    # fudge for missing UTRs

    if locus not in utr_loci:
        del v_coord_map['V-UTR']

    for coord_map in [v_coord_map, d_coord_map, j_coord_map, c_coord_map]:
        for el in coord_map.keys():
            if el not in output_headers:
                output_headers.append(el)
                output_headers.append(el + '_CIGAR')

    if locus not in utr_loci:
        output_headers.append('V-UTR')
        output_headers.append('V-UTR' + '_CIGAR')

    print("Reading assembly")
    assemblies = read_sequence_from_fasta(args.ref_seq, args.ref_seq+'.fai', refname)
    print("Reading BED files")
    beds = read_beds(refname, sense, args.bed_dir, assemblies, locus)
    rows = []

    if args.bam:
        # Process the single BAM file
        project = args.project
        subject = args.subject
        sample_name = args.sample
        sample_file = args.bam
        forced_haplotype = None
        print(f"{project} {subject} {sample_name} {sample_file}")
        rows.extend(process_sample(locus, refname, project, subject, sample_file, sample_name, beds, args.debug_gene, sense, forced_haplotype, args.list_all_contigs))
    else:
        # we expect to see project/subject/sample, with one or more bam files in the sample directory, or project/subject, with bam files in the subject directory
        projlist = glob.glob(f"{args.bam_dirs}/*")
        for proj_item in projlist:
            if os.path.isdir(proj_item):
                project = os.path.basename(proj_item)
                for subj_item in glob.glob(f"{proj_item}/*"):
                    if os.path.isdir(subj_item):
                        subject = os.path.basename(subj_item)
                        bamlist = glob.glob(f"{subj_item}/*.bam")
                        if bamlist:
                            sample_id = subject
                            process_bams(args, locus, refname, sense, beds, rows, project, subject, subj_item, sample_id, bamlist, args.list_all_contigs)
                        else:
                            for sample_item in glob.glob(f"{subj_item}/*"):
                                if os.path.isdir(sample_item):
                                    sample_id = os.path.basename(sample_item)
                                    bamlist = glob.glob(f"{sample_item}/*.bam")
                                    process_bams(args, locus, refname, sense, beds, rows, project, subject, sample_item, sample_id, bamlist, args.list_all_contigs)
                

    headers = []
    const_headers = []

    for row in rows:
        for k in row.keys():
            if 'C_' in k:
                if k not in const_headers:
                    const_headers.append(k)
            elif k not in headers:
                headers.append(k)

    const_headers.sort()
    headers.extend(const_headers)

    with open(args.outfile, 'w') as fo:
        writer = csv.DictWriter(fo, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def process_bams(args, locus, refname, sense, beds, rows, project, subject, sample_item, sample_id, bamlist, list_all_contigs):
    item_index = 1
    processed_haplotypes = {}   # store index assigned to any split haplotype files
    for sample_file in bamlist:
        forced_haplotype = None

        sn = f"{subject}_{sample_id}" if subject != sample_id else subject

        if '.1.bam' in sample_file or '.2.bam' in sample_file:
            unh = sample_file.replace('.1.bam', '.bam').replace('.2.bam', '.bam')
            forced_haplotype = 1 if '.1.bam' in sample_file else 2
            if unh not in processed_haplotypes:
                processed_haplotypes[unh] = item_index
                item_index += 1
            sample_name = f"{sn}_{processed_haplotypes[unh]}" if processed_haplotypes[unh] > 1 else f"{sn}"
        else:
            sample_name = f"{sn}_{item_index}" if item_index > 1 else f"{sn}"
            item_index += 1

        print(f"{project} {subject} {sample_name} {sample_file}")

        rows.extend(process_sample(locus, refname, project, subject, sample_file, sample_name, beds, args.debug_gene, sense, forced_haplotype, list_all_contigs))
            

if __name__ == '__main__':
    main()

