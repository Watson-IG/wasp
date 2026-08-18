import os
import argparse
import concurrent.futures
import subprocess
import pandas as pd

def _check_csv_output(csv_path, step_name, locus, loci):
    """Check that a pipeline step produced a non-empty CSV (not just headers).

    Prints a warning to stdout if the file is missing or headers-only.
    """
    label = f"{locus}/{loci}"
    if not os.path.isfile(csv_path):
        print(f"ANNOTATION_WARNING: {step_name} did not produce output file for {label}: {csv_path}")
        return
    try:
        df = pd.read_csv(csv_path)
        if len(df) == 0:
            print(f"ANNOTATION_WARNING: {step_name} produced empty table (headers only) for {label}: {csv_path}")
    except pd.errors.EmptyDataError:
        print(f"ANNOTATION_WARNING: {step_name} produced completely empty file for {label}: {csv_path}")
    except Exception as e:
        print(f"ANNOTATION_WARNING: {step_name} output could not be read for {label}: {csv_path} ({e})")

def process_locus(sample_id, bam_path, locus, loci, ref, bed_dir, allele_ref_dir, outdir):
    output_base_dir = outdir
    if loci == 'ighc':
        output_dir = f"{output_base_dir}/annotations/{sample_id}/IGHC"
    else:
        output_dir = f"{output_base_dir}/annotations/{sample_id}/{locus}"

    os.makedirs(output_dir, exist_ok=True)

    # Make gene file
    make_gene_file_output = f"{output_dir}/{sample_id}_make_gene_file.csv"
    make_command = f"/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/make_gene_file.py -l --sample {sample_id} {locus} {loci} \"+-\" {bed_dir} {ref} {make_gene_file_output} -b {bam_path}"
    subprocess.run(make_command, shell=True, check=True)
    _check_csv_output(make_gene_file_output, "make_gene_file", locus, loci)

    # Import from assembly
    import_from_assembly_output = f"{output_dir}/{sample_id}_make_gene_file_imported.csv"
    import_command = f"/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/import_from_assemblies.py {locus} {loci} \"+-\" {make_gene_file_output} {ref} {bed_dir} {allele_ref_dir} {import_from_assembly_output}"
    subprocess.run(import_command, shell=True, check=True)
    _check_csv_output(import_from_assembly_output, "import_from_assemblies", locus, loci)
    return import_from_assembly_output

def process_sample(sample_id, bam_path, ref, bed_dir, allele_ref_dir, outdir):
    loci_options = [
        ('IGH', 'igh'),
        ('IGH', 'ighc'),
        ('IGK', 'chr2'),
        ('IGL', 'chr22'),
        ('TRB', 'trb'),
        ('TRG', 'chr7'),
        ('TRD', 'chr14'),
        ('TRA', 'chr14')
    ]
    #output_base_dir = os.path.dirname(os.path.dirname(bam_path))
    #merged_outfile = f"{output_base_dir}/annotations/{sample_id}/merged_annotations.csv"

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_locus, sample_id, bam_path, locus, loci, ref, bed_dir, allele_ref_dir, outdir) for locus, loci in loci_options]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    

def main(sample_id, bam_path, ref, bed_dir, allele_ref_dir, outdir):
    process_sample(sample_id, bam_path, ref, bed_dir, allele_ref_dir, outdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process samples using make_gene_file.py and import_from_assemblies.py')
    parser.add_argument('sample_id', help='sample ID')
    parser.add_argument('bam_path', help='BAM path')
    parser.add_argument('ref', help='reference path')
    parser.add_argument('bed_dir', help='bed file dir path')
    parser.add_argument('allele_ref_dir', help='allele_ref.fasta dir path')
    parser.add_argument('outdir', help='output dir root path')
    args = parser.parse_args()
    main(args.sample_id, args.bam_path, args.ref, args.bed_dir, args.allele_ref_dir, args.outdir)
