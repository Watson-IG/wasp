#!/bin/bash
set -x

# Parse options
mode="combined"
species="Homo_sapiens"
locus_fasta_args=()
cut_distance=""
motif_dir=""

usage() {
  echo "Usage: $0 -f <config> -s <sample> -r <reads> [options]"
  echo ""
  echo "Required Options:"
  echo "  -f, --config       Path to the configuration file"
  echo "  -s, --sample       Name of the sample"
  echo "  -r, --reads        Path to the hifi reads (pacbio BAM, FASTQ, or FASTA)"
  echo ""
  echo "Optional Options:"
  echo "  -m, --mode         Analysis mode: ref, denovo, or combined (default: combined)"
  echo "  --species          Species name (default: Homo_sapiens)"
  echo "  --motif_dir        Optional path to motif directory for unsupported species"
  echo "  --locus_fasta      User-supplied per-locus assembly FASTA. Format: LOCUS=FILE"
  echo "                     (e.g. --locus_fasta IGH=/path/to/igh.fasta). Can be specified"
  echo "                     multiple times. When used, assembly (hifiasm) is skipped and"
  echo "                     the pipeline runs in denovo mode."
  echo "  --cut [int]        Trim contigs that extend beyond IG/TR loci boundaries."
  echo "                     Optionally specify a buffer distance in bp (default: 20000)."
  echo "  -h, --help         Show this help message and exit"
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--mode)
      mode="$2"
      shift 2
      ;;
    -f|--config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --species)
      species="$2"
      shift 2
      ;;
    --motif_dir)
      motif_dir="$2"
      shift 2
      ;;
    -s|--sample)
      sample="$2"
      shift 2
      ;;
    -r|--reads)
      ccs="$2"
      shift 2
      ;;
    --locus_fasta)
      locus_fasta_args+=("$2")
      shift 2
      ;;
    --cut)
      if [[ -n "$2" && ! "$2" == -* ]]; then
        cut_distance="$2"
        shift 2
      else
        cut_distance="20000"
        shift 1
      fi
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*|--*)
      echo "Unknown option $1"
      usage
      exit 1
      ;;
    *)
      echo "Unknown argument $1"
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIG_FILE" || ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Config file not found or not specified! (Use -f or --config)"
    echo ""
    usage
    exit 1
else
    source "$CONFIG_FILE"
fi

# When user supplies --locus_fasta, force denovo mode and skip assembly
user_assembly=false
if [[ ${#locus_fasta_args[@]} -gt 0 ]]; then
    user_assembly=true
    # Validate each LOCUS=FILE pair
    for lf in "${locus_fasta_args[@]}"; do
        if [[ "$lf" != *"="* ]]; then
            echo "Error: Invalid --locus_fasta format: $lf (expected LOCUS=FILE, e.g. IGH=/path/to/igh.fasta)"
            exit 1
        fi
        lf_file="${lf#*=}"
        if [[ ! -f "$lf_file" ]]; then
            echo "Error: --locus_fasta file not found: $lf_file"
            exit 1
        fi
    done
fi

if [[ -z "$sample" || -z "$ccs" ]]; then
    echo "Error: Missing required arguments! Please provide -s (sample) and -r (reads)."
    echo ""
    usage
    exit 1
fi

outdir=$PWD/run_wasp/${sample}
mkdir -p $outdir
if [ -d "$PWD/statistics" ]; then
    stats=true  
else
    stats=false
fi


if [[ "${ccs}" == *.bam ]]; then
    samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
elif [[ "${ccs}" == *.fastq || "${ccs}" == *.fq || "${ccs}" == *.fastq.gz || "${ccs}" == *.fq.gz ]]; then
    seqtk seq -a ${ccs} > ${outdir}/reads.fasta
else
    cp ${ccs} ${outdir}/reads.fasta
fi
reads="${outdir}/reads.fasta"
# Processing steps
config_base=$(basename $CONFIG_FILE)
cp $CONFIG_FILE $outdir/$config_base
echo "INPUT_CCS=$(realpath $ccs)" >> $outdir/$config_base
echo "RUN_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')" >> $outdir/$config_base
cp /opt/wasp/scripts/qc/container.yml $outdir

# CCS-to-reference coverage (always runs — only needs reads + reference)
bash /opt/wasp/scripts/qc/cov.sh "${sample}" "${ccs}" "${reference_fasta}" "${bed_dir}/IG_loci.bed" "${threads}"

motif_dir_arg=()
if [[ -n "$motif_dir" ]]; then
    motif_dir_arg=("-motif_dir" "$motif_dir")
fi

if [[ "$user_assembly" == true ]]; then
    # -----------------------------------------------------------------------
    # User-supplied assembly path: skip hifiasm, directly map and/or run digger
    # -----------------------------------------------------------------------
    echo "Running with user-supplied locus FASTAs — skipping hifiasm assembly"

    # Build --locus_fasta args for run_digger.py
    locus_fasta_py_args=("--no-blast")
    cat_fasta="${outdir}/full_asm_for_digger.fasta"
    rm -f "$cat_fasta"
    for lf in "${locus_fasta_args[@]}"; do
        locus_fasta_py_args+=("--locus_fasta" "$lf")
        cat "${lf#*=}" >> "$cat_fasta"
    done

    fofn="${outdir}/fofn.tsv"
    mkdir -p "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq"
    bam="${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam"

    if [[ "$mode" == "ref" || "$mode" == "combined" ]]; then
        echo "Mapping user-supplied contigs to reference genome..."
        minimap2 -x "${minimap_option:-asm20}" -t "${threads}" --secondary=yes -L -a "${reference_fasta}" "${cat_fasta}" > "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sam"
        samtools view -Sbh "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sam" > "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.bam"
        samtools sort -@ "${threads}" "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.bam" -o "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.tmp.sorted.bam"
        samtools index "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.tmp.sorted.bam"
        bedtools intersect -abam "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.tmp.sorted.bam" -b "${bed_dir}/IG_loci.bed" > "$bam"
        samtools index "$bam"
        cp "$cat_fasta" "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta"
        rm -f "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sam" "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.bam" "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.tmp.sorted.bam"

        bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}" "${bam}"
    fi
else
    # -----------------------------------------------------------------------
    # Standard path: hifiasm assembly + full pipeline
    # -----------------------------------------------------------------------
    bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}"
    fofn="${outdir}/fofn.tsv"
    bash /opt/wasp/scripts/hifi-mapping/pipeline.sh "${outdir}" "${ccs}" "${threads}" "${sample}" "${reference_fasta}" "${minimap_option}" "${bed_dir}" "${cut_distance}" "${allele_ref_dir}"
    bam="${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam"
fi

if [[ "$mode" == "ref" || "$mode" == "combined" ]]; then
    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${bam} ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}
fi

if [[ "$mode" == "denovo" || "$mode" == "combined" ]]; then
    if [[ "$user_assembly" == false ]]; then
        cat ${outdir}/break_at_soft_clip/1_hifi_asm.fasta ${outdir}/break_at_soft_clip/2_hifi_asm.fasta > ${outdir}/full_asm_for_digger.fasta
        locus_fasta_py_args=()
    fi
    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/run_digger.py -species "${species}" -allele_ref_dir "${allele_ref_dir}" -reads "${outdir}/reads.fasta" -minimap_option "${ccs_minimap_option}" -threads "${threads}" "${locus_fasta_py_args[@]}" "${motif_dir_arg[@]}" "${outdir}"
fi

/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/get_asm_stats.py  ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm.stats

if [[ "$mode" == "ref" || "$mode" == "combined" ]]; then
    samtools stats ${bam} > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.stats
    samtools flagstat ${bam} > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats
    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/get_read_support_VDJs.py ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed ${threads} ${outdir} ${ccs_minimap_option}
fi

# In combined mode, merge the reference-guided and digger allele tables per locus
if [[ "$mode" == "combined" ]]; then
    loci_list=("IGH" "IGK" "IGL" "TRA" "TRB" "TRD" "TRG")
    mkdir -p "${outdir}/combined_alleles"
    for locus in "${loci_list[@]}"; do
        ref_csv="${outdir}/read_support/${sample}/imported_genes/${locus}/${sample}_make_gene_file_imported_with_read_support.csv"
        digger_csv="${outdir}/digger_read_support/${locus}_digger_read_support.csv"
        combined_csv="${outdir}/combined_alleles/${sample}_${locus}_combined_alleles.csv"
        if [[ -f "$ref_csv" && -f "$digger_csv" ]]; then
            /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/merge_allele_tables.py "$ref_csv" "$digger_csv" "$combined_csv"
        elif [[ -f "$ref_csv" ]]; then
            echo "Warning: No digger results for ${locus}, skipping merge."
        elif [[ -f "$digger_csv" ]]; then
            echo "Warning: No reference-guided results for ${locus}, skipping merge."
        fi
    done
fi

if [[ "$mode" == "ref" || "$mode" == "combined" ]]; then
    bash /opt/wasp/scripts/annotation/get_vcf/final_vcf.sh ${sample} ${bam} ${reference_fasta} ${threads} ${annoconfig} ${bed_dir} "${outdir}/ccs_cov/ccs_to_ref.sorted.bam"
fi

/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/plotReadLengths.py ${outdir}/reads.fasta ${outdir}/${sample}_readLengthHistogram.png
bash /opt/wasp/scripts/qc/move_to_results.sh "${sample}" "${outdir}" "${threads}" "${config_base}" "${mode}"
/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/getFasta.py --directory $PWD/results/${sample}/alleles --mode ${mode}

# Sanity check for missing genes against gene.bed
/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/check_missing_genes.py "$PWD/results/${sample}/alleles" "${bed_dir}/gene.bed"

#if [[ $stats == true ]]; then
 #   bash /opt/wasp/scripts/qc/move_statistics_to_results.sh "${sample}" "${outdir}"
#fi