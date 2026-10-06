#!/usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser(description="Deduplicate a FASTA file by exact header and sequence combination.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_fasta", help="Output FASTA file")
    args = parser.parse_args()

    seen = set()
    with open(args.input_fasta, "r") as f_in, open(args.output_fasta, "w") as f_out:
        current_header = None
        current_seq = []
        
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    seq_str = "".join(current_seq)
                    if (current_header, seq_str) not in seen:
                        seen.add((current_header, seq_str))
                        f_out.write(f"{current_header}\n{seq_str}\n")
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header is not None:
            seq_str = "".join(current_seq)
            if (current_header, seq_str) not in seen:
                seen.add((current_header, seq_str))
                f_out.write(f"{current_header}\n{seq_str}\n")

if __name__ == "__main__":
    main()
