#!/bin/bash
set -e -x

outdir=$1
ccs=$2
threads=$3
sample=$4
reffn=$5
minimap_option=$6
bed_dir=$7

# --- NEW FUNCTION: BLAST and Split Contigs ---
function split_assembly_contigs {
    local fasta_in=$1
    local work_dir=$2
    local base_name=$(basename "${fasta_in}" .fasta)
    local blast_out="${work_dir}/${base_name}_blast_flanks.out"
    local split_fa="${work_dir}/${base_name}.split.fa"

    echo "Running BLAST/Split on ${fasta_in}..."

    # Canonical flanks path
    # Using outfmt "6 std qlen" for coverage filtering
    blastn -query /opt/wasp/scripts/refs/5-10_left_right_flanks.fasta \
           -subject "${fasta_in}" \
           -out "${blast_out}" \
           -outfmt "6 std qlen"

    # Python script to split contigs
    python3 - "${blast_out}" "${fasta_in}" "${split_fa}" << 'PY'
import sys
import pysam
import os

# --- Settings ---
MIN_IDENTITY = 99.0
MIN_COVERAGE = 0.90  # 90% of the flank length must align

blast_file = sys.argv[1]
contigs_file = sys.argv[2]
split_file = sys.argv[3]

def get_side(query_name):
    q = query_name.lower()
    return "left" if "left" in q else ("right" if "right" in q else None)

# --- Parse BLAST (Expects outfmt "6 std qlen") ---
best_hits = {"left": {}, "right": {}}

if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
    with open(blast_file) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 13: continue 
            
            qname, sname = cols[0], cols[1]
            pident = float(cols[2])
            qstart, qend = int(cols[6]), int(cols[7])
            sstart, send = int(cols[8]), int(cols[9])
            bitscore = float(cols[11])
            qlen = int(cols[12])

            side = get_side(qname)
            
            if not side or pident < MIN_IDENTITY:
                continue

            alignment_len = abs(qend - qstart) + 1
            if (alignment_len / qlen) < MIN_COVERAGE:
                continue

            current = best_hits[side].get(sname)
            if current is None or bitscore > current["bits"]:
                best_hits[side][sname] = {
                    "sstart": sstart, "send": send, "bits": bitscore
                }

# --- Process Sequences ---
fa = pysam.FastaFile(contigs_file)
contig_lens = {r: fa.get_reference_length(r) for r in fa.references}
targets = set(best_hits["left"]).intersection(best_hits["right"])
changed = False

with open(split_file, "w") as out:
    for ref in fa.references:
        seq = fa.fetch(ref)
        length = contig_lens[ref]

        if ref not in targets:
            out.write(f">{ref}\n{seq}\n")
            continue

        L_pos = best_hits["left"][ref]["sstart"]
        R_pos = best_hits["right"][ref]["send"]

        start_cut, end_cut = sorted([L_pos, R_pos])
        start_cut = max(1, start_cut)
        end_cut = min(length, end_cut)

        # Part 1: Upstream
        if start_cut > 1:
            out.write(f">{ref}__1_{start_cut-1}\n{seq[:start_cut-1]}\n")
        
        # Part 2: Target Region
        out.write(f">{ref}__{start_cut}_{end_cut}\n{seq[start_cut-1:end_cut]}\n")
        
        # Part 3: Downstream
        if end_cut < length:
            out.write(f">{ref}__{end_cut+1}_{length}\n{seq[end_cut:]}\n")
            
        changed = True

print("SPLIT" if changed else "NO_SPLIT")
PY

    # If split file exists and is not empty, replace original
    # We check if python printed "SPLIT" or just look for file existence logic
    if [ -s "${split_fa}" ]; then
        echo "Applying splits to ${fasta_in}"
        mv "${split_fa}" "${fasta_in}"
    else
        echo "No splits required for ${fasta_in}"
        rm -f "${split_fa}" 2>/dev/null || true
    fi
}

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
    
    # Removed BLAST/Split/Re-align logic here as it is now done upstream on Hifiasm output
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

        # --- UPDATED: Split contigs BEFORE indexing and downstream processing ---
        split_assembly_contigs "${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta" "${outdir}/hifiasm"
        # ------------------------------------------------------------------------

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