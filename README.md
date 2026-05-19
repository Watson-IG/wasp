# WASP README
![f706849b-0f9b-47fe-b0da-a0d612153eaa-1](https://github.com/user-attachments/assets/1813908b-760e-432c-b6e6-c0e9bd26cd2a)

## Overview

This toolset facilitates the alignment and processing of immune receptor genomics data using a containerized pipeline. The pipeline incorporates user-provided allele reference data and outputs structured results for downstream analysis.

## Prerequisites

1. **Allele Directory** (`allele_dir`):  
   - A directory containing reference FASTA files required for importing closest allele matches.  
   - You can use any appropriate allele references. The following resources are suggestions:  
     - [Reference FASTA](http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta)  
     - [Reference FASTA (tar.gz)](http://immunogenomics.louisville.edu/wasp/ref.tar.gz)  
   - Details on a sample reference setup can be found in this [GitHub Repository](https://github.com/Watson-IG/immune_receptor_genomics/tree/main).

2. **Container**:  
   - The container image required for running the pipeline: [WASP Container (SIF)](http://immunogenomics.louisville.edu/wasp/wasp-241023.sif).

## Usage

Run the pipeline using Singularity with the following command:

    singularity exec --bind /home:/home ${container_path} /opt/wasp/scripts/wasp_align-pacbio.sh ${CONFIG_FILE} ${sample} ${ccs}

## Output

1. **Intermediary Files**:  
   These are stored in `${PWD}/run_wasp/${sampleId}` and include outputs from steps such as hifiasm and pre-soft-clipping.

2. **Final Results**:  
   The final results are stored in `${PWD}/results/${sampleId}`. Key directories include:
   - `alignments/`: Sorted BAM files and indices.
   - `alleles/`: Annotated allele files with read support.
   - `reads/`: Filtered contigs and CCS reads.
   - `stats/`: Assembly statistics and depth information.
   - `variants/`: Annotated VCF files.

## Citation

If this repository contributes to your work, particularly in the curation of immune loci genes and alleles (IG/TR) or in producing read-support for assembly validation or curated genes, please cite the following work:

https://doi.org/10.5281/zenodo.20141874
