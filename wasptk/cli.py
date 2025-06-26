import argparse
from .read_support import compute_read_support


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute read support from a BAM")
    parser.add_argument("allele_table", help="Path to allele_annotation.csv")
    parser.add_argument("bam", help="Mapped reads BAM")
    parser.add_argument("output", help="Output CSV")
    parser.add_argument("--contig-col", default="contig", help="Column for contig name")
    parser.add_argument("--start-col", default="start", help="Column for start position")
    parser.add_argument("--end-col", default="end", help="Column for end position")
    parser.add_argument("--gene-col", default="gene", help="Column for gene name")

    args = parser.parse_args()
    compute_read_support(
        args.allele_table,
        args.bam,
        args.output,
        contig_col=args.contig_col,
        start_col=args.start_col,
        end_col=args.end_col,
        gene_col=args.gene_col,
    )


if __name__ == "__main__":
    main()
