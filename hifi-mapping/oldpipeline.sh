#!/bin/bash
set -e -x

outdir=$1
ccs=$2
threads=$3
sample=$4
reffn=$5
minimap_option=$6
bed_dir=$7

function align_with_minimap2 {
    fasta=$1
    prefix=$2
    reffn=$3
    threads=$4
    minimap2 \
	-t ${threads} --secondary=yes -L -a ${reffn} \
	${fasta} > ${prefix}.sam    
    samtools view -Sbh ${prefix}.sam > ${prefix}.bam   
    samtools sort -@ ${threads} ${prefix}.bam -o ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    rm -f ${prefix}.sam
    rm -f ${prefix}.bam
}

function align_with_minimap2_asm20 {
    fasta=$1
    prefix=$2
    reffn=$3
    threads=$4
    minimap_option=$5
    minimap2 -x ${minimap_option} \
	-t ${threads} --secondary=yes -L -a ${reffn} \
	${fasta} > ${prefix}.sam    
    samtools view -Sbh ${prefix}.sam > ${prefix}.bam   
    samtools sort -@ ${threads} ${prefix}.bam -o ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    rm -f ${prefix}.sam
    rm -f ${prefix}.bam
}

function align_and_process {
    sample=$1
    dir=$2
    outdir=${dir}/merged_bam/final_asm20_to_ref_with_secondarySeq
    mkdir -p $outdir
    minimap2 -x asm20 --secondary=yes -t ${threads} -L -a ${reffn} ${dir}/merged_bam/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ ${threads} ${outdir}/${sample}.bam -o ${outdir}/${sample}.tmp.sorted.bam
    samtools index ${outdir}/${sample}.tmp.sorted.bam
    bedtools intersect -abam ${outdir}/${sample}.tmp.sorted.bam -b ${bed_dir}/IG_loci.bed > ${outdir}/${sample}.sorted.bam
    samtools index ${outdir}/${sample}.sorted.bam
    samtools view ${outdir}/${sample}.sorted.bam | awk '{ print ">"$1"\n"$10 }' > ${outdir}/contigs.fasta

}
function merge_and_rmdup {
    sample=$1
    dir=$2
    outdir=${dir}/merged_bam/final_asm20_to_ref_with_secondarySeq
    mkdir -p $outdir

    # Add read group information using temporary files
    samtools addreplacerg -r "ID:1 HPL:hap1" -o ${dir}/break_at_soft_clip/temp1.bam ${dir}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam
    mv ${dir}/break_at_soft_clip/temp1.bam ${dir}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam
    samtools addreplacerg -r "ID:2 HPL:hap2" -o ${dir}/break_at_soft_clip/temp2.bam ${dir}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam
    mv ${dir}/break_at_soft_clip/temp2.bam ${dir}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam
    
    # Merge BAM files with read group tags
    samtools merge -f ${outdir}/merged.bam ${dir}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam ${dir}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam
    
    # Sort and index the merged file
    samtools sort -@ ${threads} ${outdir}/merged.bam -o ${outdir}/${sample}.tmp.sorted.bam
    samtools index ${outdir}/${sample}.tmp.sorted.bam
    
    # Remove duplicates and index the final BAM file
    samtools markdup -r ${outdir}/${sample}.tmp.sorted.bam ${outdir}/${sample}.rmdup.bam
    samtools index ${outdir}/${sample}.rmdup.bam
    
    bedtools intersect -abam ${outdir}/${sample}.rmdup.bam -b ${bed_dir}/IG_loci.bed > ${outdir}/${sample}.sorted.bam
    samtools index ${outdir}/${sample}.sorted.bam
    samtools view ${outdir}/${sample}.sorted.bam | awk '{ print ">"$1"\n"$10 }' > ${outdir}/contigs.fasta
}




if [ ! -s ${outdir}/reads.fasta.fai ]
then
   # samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
    samtools faidx ${outdir}/reads.fasta
fi

if [ ! -s ${outdir}/hifiasm/asm.bp.hap2.p_ctg.fasta.fai ]
then
    

    if [ ! -f ${outdir}/hifiasm/asm.bp.hap2.p_ctg.gfa ]
    then
		mkdir -p ${outdir}/hifiasm
        hifiasm \
        -o ${outdir}/hifiasm/asm \
        -t ${threads} \
        ${outdir}/reads.fasta
    fi

    for i in p #r
    do
        gfatools \
            gfa2fa \
            ${outdir}/hifiasm/asm.bp.${i}_utg.gfa > \
            ${outdir}/hifiasm/asm.bp.${i}_utg.fasta
    done

    for i in 1 2
    do
        gfatools \
            gfa2fa \
            ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.gfa > \
            ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta

        # Add suffix "_hap${i}" to each header in the fasta file
        awk '/^>/{print $0 "_hap'${i}'"} !/^>/{print}' ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta > ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.modified.fasta

        # Move modified file to original name if needed
        mv ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.modified.fasta ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta

        samtools faidx ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta
    done

fi


for i in 1 2
do
    fn=asm.bp.hap${i}.p_ctg
    if [ ! -s ${outdir}/hifiasm/${fn}_to_ref.sorted.bam.bai ]
    then
	align_with_minimap2_asm20 \
	    ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta \
	    ${outdir}/hifiasm/${fn}_to_ref \
	    ${reffn} ${threads} "${minimap_option}"
    fi
done

for i in p #r
do
    fn=asm.bp.${i}_utg
    if [ ! -s ${outdir}/hifiasm/${fn}_to_ref.sorted.bam.bai ]
    then
	align_with_minimap2_asm20 \
	    ${outdir}/hifiasm/${fn}.fasta \
	    ${outdir}/hifiasm/${fn}_to_ref \
	    ${reffn} \
	    ${threads} "${minimap_option}"
    fi
done


mkdir -p ${outdir}/break_at_soft_clip

for i in 1 2
do
    bam=${outdir}/hifiasm/asm.bp.hap${i}.p_ctg_to_ref.sorted.bam

    if [ ! -s ${outdir}/break_at_soft_clip/${i}_hifi_asm_to_ref.sorted.bam ]
    then
        python /opt/wasp/scripts/hifi-mapping/extract_soft_clip_seq.py \
        ${bam} > ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta

        samtools faidx ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta

       # align_with_minimap2 \
        #${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta \
        #${outdir}/break_at_soft_clip/${i}_hifi_asm_to_ref \
        #${reffn} \
        #${threads}
    fi

    if [ ! -s ${outdir}/break_at_soft_clip/${i}_asm20_hifi_asm_to_ref.sorted.bam ]
    then
        align_with_minimap2_asm20 \
        ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta \
        ${outdir}/break_at_soft_clip/${i}_asm20_hifi_asm_to_ref \
        ${reffn} \
        ${threads} "${minimap_option}"
    fi
done
merge_and_rmdup $sample $outdir
#align_and_process $sample $outdir

