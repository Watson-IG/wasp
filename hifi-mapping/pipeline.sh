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
    blastn -query /opt/wasp/scripts/refs/5-10_left_right_flanks.fasta -subject ${outdir}/contigs.fasta  -out ${outdir}/contigs_blast_flanks.out -outfmt 7
        # split contigs if both flanks hit at >=99% identity, then realign contigs and regenerate ${sample}.sorted.bam
    python3 - <<'PY'
import os, pysam, sys
blast = os.path.join(os.environ["outdir"], "contigs_blast_flanks.out")
contigs = os.path.join(os.environ["outdir"], "contigs.fasta")
splitfa = os.path.join(os.environ["outdir"], "contigs.split.fa")

def side(q):
    q=q.lower()
    return "left" if "left" in q else ("right" if "right" in q else None)

best={"left":{}, "right":{}}
with open(blast) as fh:
    for line in fh:
        if not line or line.startswith("#"): continue
        p=line.rstrip("\n").split("\t")
        if len(p)<12: continue
        s=side(p[0]); 
        if not s or float(p[2])<99.0: continue
        sseq=p[1]; sstart=int(p[8]); send=int(p[9]); bits=float(p[11])
        cur=best[s].get(sseq)
        if cur is None or bits>cur["bits"]:
            best[s][sseq]={"sstart":sstart,"send":send,"bits":bits}

fa=pysam.FastaFile(contigs)
lens={r:fa.get_reference_length(r) for r in fa.references}
pairs=set(best["left"]).intersection(best["right"])
changed=False
with open(splitfa,"w") as out:
    def w(name, seq):
        for i in range(0,len(seq),60):
            out.write(seq[i:i+60]+"\n")
    for r in fa.references:
        if r not in pairs:
            out.write(f">{r}\n"); w(r, fa.fetch(r)); continue
        L=best["left"][r]["sstart"]; R=best["right"][r]["send"]
        if L>R: L,R=R,L
        L=max(1,L); R=min(lens[r],R)
        out.write(f">{r}__1_{L-1}\n"); 
        if L>1: w(r, fa.fetch(r, 0, L-1))
        else: out.write("\n")
        out.write(f">{r}__{L}_{R}\n"); w(r, fa.fetch(r, L-1, R))
        out.write(f">{r}__{R+1}_{lens[r]}\n");
        if R<lens[r]: w(r, fa.fetch(r, R, lens[r]))
        else: out.write("\n")
        changed=True
print("SPLIT" if changed else "NO_SPLIT")
PY

    if [ "$(tail -n1 ${outdir}/contigs.split.fa 2>/dev/null)" != "" ]; then
        mv ${outdir}/contigs.split.fa ${outdir}/contigs.fasta
        samtools faidx ${outdir}/contigs.fasta
    else
        rm -f ${outdir}/contigs.split.fa 2>/dev/null || true
    fi

    # realign contigs to reference using the same settings and regenerate outputs
    minimap2 -x asm20 --secondary=yes -t ${threads} -L -a ${reffn} ${outdir}/contigs.fasta > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ ${threads} ${outdir}/${sample}.bam -o ${outdir}/${sample}.tmp.sorted.bam
    samtools index ${outdir}/${sample}.tmp.sorted.bam
    bedtools intersect -abam ${outdir}/${sample}.tmp.sorted.bam -b ${bed_dir}/IG_loci.bed > ${outdir}/${sample}.sorted.bam
    samtools index ${outdir}/${sample}.sorted.bam

}

function merge_and_rmdup {
    local sample=$1
    local dir=$2
    local outdir=${dir}/merged_bam/final_asm20_to_ref_with_secondarySeq
    mkdir -p "$outdir"

    # Add valid read-group tags
    samtools addreplacerg -r "ID:1" -r "SM:hap1" -o "${dir}/break_at_soft_clip/temp1.bam" "${dir}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam"
    mv "${dir}/break_at_soft_clip/temp1.bam" "${dir}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam"

    samtools addreplacerg -r "ID:2" -r "SM:hap2" -o "${dir}/break_at_soft_clip/temp2.bam" "${dir}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam"
    mv "${dir}/break_at_soft_clip/temp2.bam" "${dir}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam"
    
    # Merge, sort, index
    samtools merge -f "${outdir}/merged.bam" \
        "${dir}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam" \
        "${dir}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam"
    samtools sort -@ "${threads}" "${outdir}/merged.bam" -o "${outdir}/${sample}.tmp.sorted.bam"
    samtools index "${outdir}/${sample}.tmp.sorted.bam"
    
    # Duplicate flagging
    in_bam="${outdir}/${sample}.tmp.sorted.bam"
    out_bam="${outdir}/${sample}.rmdup.bam"
    python3 - "$in_bam" "$out_bam" << 'PYCODE'
import sys, pysam
inBam, outBam = sys.argv[1], sys.argv[2]
with pysam.AlignmentFile(inBam, "rb") as bamIn, \
     pysam.AlignmentFile(outBam, "wb", template=bamIn) as bamOut:
    readCache = {}
    currentRef, lastPos = None, -1
    for read in bamIn.fetch(until_eof=True):
        if read.reference_id != currentRef or read.reference_start != lastPos:
            readCache.clear()
            currentRef, lastPos = read.reference_id, read.reference_start
        key = (read.reference_id, read.reference_start, read.is_reverse,
               read.query_sequence, read.cigarstring)
        read.is_duplicate = key in readCache
        if not read.is_duplicate:
            readCache[key] = True
        bamOut.write(read)
PYCODE
    samtools index "${out_bam}"
    
    # Filter to IG loci, fasta export
    bedtools intersect -abam "${out_bam}" -b "${bed_dir}/IG_loci.bed" > "${outdir}/${sample}.sorted.bam"
    samtools index "${outdir}/${sample}.sorted.bam"
    samtools view "${outdir}/${sample}.sorted.bam" | awk '{ print ">"$1"\n"$10 }' > "${outdir}/contigs.fasta"

    # Canonical flanks path
    blastn -query /opt/wasp/scripts/refs/5-10_left_right_flanks.fasta \
           -subject "${outdir}/contigs.fasta" \
           -out "${outdir}/contigs_blast_flanks.out" -outfmt 7

    # Split contigs if both flanks hit at >=99% identity
    python3 - "${outdir}/contigs_blast_flanks.out" "${outdir}/contigs.fasta" "${outdir}/contigs.split.fa" << 'PY'
import sys, pysam
blast, contigs, splitfa = sys.argv[1], sys.argv[2], sys.argv[3]

def side(q):
    q=q.lower()
    return "left" if "left" in q else ("right" if "right" in q else None)

best={"left":{}, "right":{}}
with open(blast) as fh:
    for line in fh:
        if not line or line.startswith("#"): continue
        p=line.rstrip("\n").split("\t")
        if len(p)<12: continue
        s=side(p[0])
        if not s or float(p[2])<99.0: continue
        sseq=p[1]; sstart=int(p[8]); send=int(p[9]); bits=float(p[11])
        cur=best[s].get(sseq)
        if cur is None or bits>cur["bits"]:
            best[s][sseq]={"sstart":sstart,"send":send,"bits":bits}

fa=pysam.FastaFile(contigs)
lens={r:fa.get_reference_length(r) for r in fa.references}
pairs=set(best["left"]).intersection(best["right"])
changed=False
with open(splitfa,"w") as out:
    def w(seq):
        for i in range(0,len(seq),60):
            out.write(seq[i:i+60]+"\n")
    for r in fa.references:
        if r not in pairs:
            out.write(f">{r}\n"); w(fa.fetch(r)); continue
        L=best["left"][r]["sstart"]; R=best["right"][r]["send"]
        if L>R: L,R=R,L
        L=max(1,L); R=min(lens[r],R)
        out.write(f">{r}__1_{L-1}\n"); w(fa.fetch(r, 0, L-1) if L>1 else "")
        out.write(f">{r}__{L}_{R}\n"); w(fa.fetch(r, L-1, R))
        out.write(f">{r}__{R+1}_{lens[r]}\n"); w(fa.fetch(r, R, lens[r]) if R<lens[r] else "")
        changed=True
print("SPLIT" if changed else "NO_SPLIT")
PY

    if [ -s "${outdir}/contigs.split.fa" ]; then
        mv "${outdir}/contigs.split.fa" "${outdir}/contigs.fasta"
        samtools faidx "${outdir}/contigs.fasta"
    else
        rm -f "${outdir}/contigs.split.fa" 2>/dev/null || true
    fi

    # Realign post-split and regenerate outputs
    minimap2 -x asm20 --secondary=yes -t "${threads}" -L -a "${reffn}" "${outdir}/contigs.fasta" > "${outdir}/${sample}.sam"
    samtools view -Sbh "${outdir}/${sample}.sam" > "${outdir}/${sample}.bam"
    samtools sort -@ "${threads}" "${outdir}/${sample}.bam" -o "${outdir}/${sample}.tmp.sorted.bam"
    samtools index "${outdir}/${sample}.tmp.sorted.bam"
    bedtools intersect -abam "${outdir}/${sample}.tmp.sorted.bam" -b "${bed_dir}/IG_loci.bed" > "${outdir}/${sample}.sorted.bam"
    samtools index "${outdir}/${sample}.sorted.bam"
}


if [ ! -s ${outdir}/reads.fasta.fai ]
then
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

        awk '/^>/{print $0 "_hap'${i}'"} !/^>/{print}' ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta > ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.modified.fasta

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

    if [ ! -s ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta ]
    then
        python /opt/wasp/scripts/hifi-mapping/extract_soft_clip_seq.py \
            ${bam} > ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta

        samtools faidx ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta
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
# align_and_process $sample $outdir
