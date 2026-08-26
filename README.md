# WASP: Immune Receptor Genomics Pipeline

![f706849b-0f9b-47fe-b0da-a0d612153eaa-1](https://github.com/user-attachments/assets/1813908b-760e-432c-b6e6-c0e9bd26cd2a)

## Overview

WASP is a containerized pipeline designed to process, assemble, and annotate PacBio HiFi reads for immune receptor loci (IG/TR). It performs read mapping, assembly (via `hifiasm`), and highly specific allele annotation.

## Prerequisites

1. **Allele Directory** (`allele_dir`):  
   - A directory containing reference FASTA files required for importing closest allele matches.  
   - Suggestions include:  
     - [Reference FASTA](http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta)  
     - [Reference FASTA (tar.gz)](http://immunogenomics.louisville.edu/wasp/ref.tar.gz)  
   - Details on a sample reference setup can be found in this [GitHub Repository](https://github.com/Watson-IG/immune_receptor_genomics/tree/main).

2. **Container**:  
   - A WASP Singularity container image (SIF) is required to run the pipeline. 

## Usage

The primary execution script is `wasp.sh`. You run it through Singularity:

```bash
singularity exec --bind /home:/home ${container_path} /opt/wasp/scripts/wasp.sh -f ${CONFIG_FILE} -s ${sample} -r ${ccs} [options]
```

### Required Arguments
- `-f, --config` : Path to the configuration file.
- `-s, --sample` : Name of the sample.
- `-r, --reads`  : Path to the HiFi reads (PacBio BAM, FASTQ, or FASTA).

### Analysis Modes (`-m, --mode`)
The pipeline supports three annotation modes (default is `combined`):
- **`ref` (Reference-guided)**: Maps the assembled contigs against known allele references to annotate genes and alleles.
- **`denovo` (De novo)**: Uses `digger` to perform *de novo* motif-based annotation directly on the assembled contigs, which is useful for discovering novel alleles or structural variations.
- **`combined`**: Runs both `ref` and `denovo` modes. The results are merged into a unified allele table, keeping reference annotations when overlapping with digger annotations.

### Bypassing Assembly with Pre-binned Contigs
If you already have assembled contigs and want to skip the built-in `hifiasm` assembly step, you can pass them directly to the pipeline using the `--locus_fasta` flag.
- **Format**: `LOCUS=FILE`
- **Example**: `--locus_fasta IGH=/path/to/igh.fasta --locus_fasta IGK=/path/to/igk.fasta`

*Note: Supplying pre-binned contigs will automatically force the pipeline into `denovo` mode, as it skips the standard assembly and reference-mapping steps, running `digger` directly on your supplied FASTA files.*

### Additional Options
- `--species` : Species name (default: `Homo_sapiens`).
- `--motif_dir` : Optional path to a custom digger motif directory (useful for unsupported species).
- `--cut [int]` : Trim contigs that extend beyond IG/TR loci boundaries. Optionally specify a buffer distance in bp (default: 20000).

## Output

1. **Intermediary Files**:  
   Stored in `${PWD}/run_wasp/${sample}` and include raw outputs from tools like `hifiasm` and `digger`.

2. **Final Results**:  
   The final results are stored in `${PWD}/results/${sample}`. Key directories include:
   - `alignments/`: Sorted BAM files and indices.
   - `alleles/`: Annotated allele tables (e.g., `_combined_annotated-alleles.csv`) with read support.
   - `reads/`: Filtered contigs and CCS reads.
   - `stats/`: Assembly statistics and depth information.
   - `variants/`: Annotated VCF files.

## Citation

If this repository contributes to your work, particularly in the curation of immune loci genes and alleles (IG/TR) or in producing read-support for assembly validation or curated genes, please cite the following work:

https://doi.org/10.5281/zenodo.20141874
