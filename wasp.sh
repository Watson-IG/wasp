#!/bin/bash
set -x

CONFIG_FILE=$1
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "Config file not found!"
    exit 1
fi

sample=$2
ccs=$3
asm=$4
outdir=${5:-$PWD/run_wasp/${sample}}

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

cp $CONFIG_FILE $outdir
cp /opt/wasp/scripts/qc/container.yml $outdir
bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}" "${asm}"
fofn="${outdir}/fofn.tsv"
bash /opt/wasp/scripts/qc/cov.sh "${sample}" "${ccs}" "${reference_fasta}" "${bed_dir}/IG_loci.bed" "${threads}"
echo "/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${asm} ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}"
/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${asm} ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}
/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/get_asm_stats.py ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm.stats
samtools stats ${asm} > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.stats
samtools flagstat ${asm} > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats
/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/get_read_support_VDJs.py ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed ${threads} ${outdir} ${ccs_minimap_option}
bash /opt/wasp/scripts/annotation/get_vcf/final_vcf.sh ${sample} ${asm} ${reference_fasta} ${threads} ${annoconfig} ${bed_dir}
#bash /opt/wasp/scripts/qc/perscov.sh "${sample}" "${outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam" "${asm}" "${bed_dir}/IG_loci.bed" "${outdir}"
/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/plotReadLengths.py ${outdir}/reads.fasta ${outdir}/${sample}_readLengthHistogram.png
bash /opt/wasp/scripts/qc/move_to_results.sh "${sample}" "${outdir}"
#if [[ $stats == true ]]; then
#    bash /opt/wasp/scripts/qc/move_statistics_to_results.sh "${sample}" "${outdir}"
#fi
