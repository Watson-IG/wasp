#!/bin/bash
set -x

# Parse options
mode="ref"
species="human"

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
    -S|--species)
      species="$2"
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
    -c|--container)
      container="$2"
      shift 2
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      echo "Unknown argument $1"
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIG_FILE" || ! -f "$CONFIG_FILE" ]]; then
    echo "Config file not found or not specified! (Use -f or --config)"
    exit 1
else
    source "$CONFIG_FILE"
fi

if [[ -z "$sample" || -z "$ccs" || -z "$container" ]]; then
    echo "Missing required arguments! Please provide -s (sample), -r (reads), and -c (container)."
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
else
    cp ${ccs} ${outdir}/reads.fasta
fi
reads="${outdir}/reads.fasta"
# Processing steps
#singularity exec ${container}
config_base=$(basename $CONFIG_FILE)
cp $CONFIG_FILE $outdir/$config_base
echo "INPUT_CCS=$(realpath $ccs)" >> $outdir/$config_base
echo "RUN_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')" >> $outdir/$config_base
cp /opt/wasp/scripts/qc/container.yml $outdir
bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}"
fofn="${outdir}/fofn.tsv"
bash /opt/wasp/scripts/qc/cov.sh "${sample}" "${ccs}" "${reference_fasta}" "${bed_dir}/IG_loci.bed" "${threads}"
bash /opt/wasp/scripts/hifi-mapping/pipeline.sh "${outdir}" "${ccs}" "${threads}" "${sample}" "${reference_fasta}" "${minimap_option}" "${bed_dir}"
if [[ "$mode" == "ref" || "$mode" == "combined" ]]; then
    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}
fi

if [[ "$mode" == "denovo" || "$mode" == "combined" ]]; then
    cat ${outdir}/break_at_soft_clip/1_hifi_asm.fasta ${outdir}/break_at_soft_clip/2_hifi_asm.fasta > ${outdir}/full_asm_for_digger.fasta
    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/run_digger.py -species "${species}" -allele_ref_dir "${allele_ref_dir}" -reads "${outdir}/reads.fasta" -minimap_option "${ccs_minimap_option}" -threads "${threads}" "${outdir}"
fi
/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/get_asm_stats.py  ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm.stats
samtools stats ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.stats
samtools flagstat ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats
if [[ "$mode" == "ref" || "$mode" == "combined" ]]; then
    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/get_read_support_VDJs.py ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed ${threads} ${outdir} ${ccs_minimap_option}
fi
bash /opt/wasp/scripts/annotation/get_vcf/final_vcf.sh ${sample} ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${reference_fasta} ${threads} ${annoconfig} ${bed_dir} "${outdir}/ccs_cov/ccs_to_ref.sorted.bam"
#bash /opt/wasp/scripts/qc/perscov.sh "${sample}" "${outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam" "${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam" "${bed_dir}/IG_loci.bed" "${outdir}"
/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/plotReadLengths.py ${outdir}/reads.fasta ${outdir}/${sample}_readLengthHistogram.png
bash /opt/wasp/scripts/qc/move_to_results.sh "${sample}" "${outdir}" "${threads}" "${config_base}"
/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/getFasta.py --directory $PWD/results/${sample}/alleles
#if [[ $stats == true ]]; then
 #   bash /opt/wasp/scripts/qc/move_statistics_to_results.sh "${sample}" "${outdir}"
#fi