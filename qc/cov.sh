#!/bin/bash
set -e -x

scratch=$PWD
sample=$1
ccs=$2
reffn=$3
#refbed=/home/zmvanw01/test-beds/sorted_region.bed
#refbed=/home/zmvanw01/projects/t_Parks/parks.bed
#refbed=/home/zmvanw01/240520-coverage.bed
refbed=$4
threads=$5

function run_get_ccs_cov {
    mkdir -p ${scratch}/run_wasp/${sample}/ccs_cov
    outd=${scratch}/run_wasp/${sample}/ccs_cov
    bam_path=${scratch}/run_wasp/${sample}/ccs_cov/ccs_to_ref.sorted.bam
    pers_bam_path=${scratch}/run_wasp/${sample}/ccs_cov/ccs_to_pers-ref.sorted.bam
    #bamtocov --regions ${refbed} --report ${outd}/${sample}_stats.tsv ${bam_path} > ${outd}/${sample}_cov-cov.bed
    mosdepth -b ${refbed} -t ${threads} ${outd}/${sample} ${bam_path}
}
function map_ccs_to_ref {
    #dir=$scratch/ccs_cov
    readspath=${scratch}/run_wasp/${sample}/reads.fasta
    outdir=${scratch}/run_wasp/${sample}/ccs_cov
    mkdir -p $outdir
    #samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
    samtools faidx ${readspath}
    minimap2 -x map-hifi --secondary=no -t "${threads}" -L -a ${reffn} ${readspath} > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ "${threads}" ${outdir}/${sample}.bam -o ${outdir}/ccs_to_ref.sorted.bam
    samtools index ${outdir}/ccs_to_ref.sorted.bam
    samtools flagstats ${outdir}/ccs_to_ref.sorted.bam > ${outdir}/ccs_to_ref-full.flagstats.txt
    samtools view -L ${refbed} -b ${outdir}/ccs_to_ref.sorted.bam > ${outdir}/ig-filtered_ccs_to_ref.sorted.bam
    samtools index ${outdir}/ig-filtered_ccs_to_ref.sorted.bam
    samtools flagstats ${outdir}/ig-filtered_ccs_to_ref.sorted.bam > ${outdir}/ig-filtered_ccs_to_ref.flagstats.txt
    rm ${outdir}/ig-filtered_ccs_to_ref.sorted.bam

    on_target=$(awk 'NR==1 {print $1}' "${outdir}/ig-filtered_ccs_to_ref.flagstats.txt")
    full=$(awk 'NR==1 {print $1}' "${outdir}/ccs_to_ref-full.flagstats.txt")
    percent_on=$(awk -v a="$on_target" -v b="$full" 'BEGIN { printf "%.2f", (100 * a / b) }')
    echo "$percent_on" > "${outdir}/percent_on_target.txt"

}
function map_ccs_to_pers_ref {
    dir=$scratch/ccs_cov
    outdir=${scratch}/run_wasp/${sample}/ccs_cov
    mkdir -p $outdir
    minimap2 -x map-hifi --secondary=no -t 12 -L -a ${reffn} $PWD/run_wasp/${sample}/read_support/${sample}/ccs_to_pers/pers_ref.fasta > ${outdir}/${sample}.pers.sam
    samtools view -Sbh ${outdir}/${sample}.pers.sam > ${outdir}/${sample}.pers.bam
    samtools sort -@ 12 ${outdir}/${sample}.pers.bam -o ${outdir}/ccs_to_pers-ref.sorted.bam
    samtools index ${outdir}/ccs_to_pers-ref.sorted.bam
}
#map_ccs_to_pers_ref
map_ccs_to_ref
run_get_ccs_cov
