#!/bin/bash

sample=$1
orig_outdir=$2
threads=$3
config_file=$4
mode=${5:-ref}
outdir=$PWD/results/${sample}

# Helper: safely move a file and create a symlink back to the original location.
# If the source file doesn't exist, print a warning and continue instead of crashing.
safe_mv_link() {
    local src="$1"
    local dst="$2"
    if [ -f "$src" ]; then
        mv "$src" "$dst"
        ln -s "$dst" "$src"
    else
        echo "MOVE_WARNING: Expected file not found, skipping: $src"
    fi
}

mkdir -p ${outdir}
mkdir -p ${outdir}/reads ${outdir}/alignments ${outdir}/variants ${outdir}/alleles ${outdir}/stats ${outdir}/metadata

# Moving metadata and config
if [ -f "${orig_outdir}/container.yml" ]; then
    mv ${orig_outdir}/container.yml ${outdir}/metadata/container.yml
    ln -s ${outdir}/metadata/container.yml ${orig_outdir}/container.yml
fi

if [ ! -z "${config_file}" ] && [ -f "${orig_outdir}/${config_file}" ]; then
    mv "${orig_outdir}/${config_file}" "${outdir}/metadata/${config_file}"
    ln -s "${outdir}/metadata/${config_file}" "${orig_outdir}/${config_file}"
fi

# Moving and creating symlinks
safe_mv_link "${orig_outdir}/reads.fasta" "${outdir}/reads/ccs-reads.fasta"

safe_mv_link "${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta" "${outdir}/reads/hifiasm_ig-filtered_contigs.fasta"

safe_mv_link "${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam" "${outdir}/alignments/${sample}_contigs-to-ref.sorted.bam"

safe_mv_link "${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam.bai" "${outdir}/alignments/${sample}_contigs-to-ref.sorted.bam.bai"

safe_mv_link "${orig_outdir}/ccs_cov/ccs_to_ref.sorted.bam" "${outdir}/alignments/${sample}_ccs-to-ref.sorted.bam"

safe_mv_link "${orig_outdir}/ccs_cov/ccs_to_ref.sorted.bam.bai" "${outdir}/alignments/${sample}_ccs-to-ref.sorted.bam.bai"

safe_mv_link "${orig_outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam" "${outdir}/alignments/${sample}_ccs-to-personal-reference.sorted.bam"

safe_mv_link "${orig_outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam.bai" "${outdir}/alignments/${sample}_ccs-to-personal-reference.sorted.bam.bai"

# Loci list
loci_list=("IGH" "IGHC" "IGK" "IGL" "TRA" "TRB" "TRD" "TRG")

for loci in "${loci_list[@]}"; do
    src="${orig_outdir}/read_support/${sample}/imported_genes/${loci}/${sample}_make_gene_file_imported_with_read_support.csv"
    dst="${outdir}/alleles/${sample}_${loci}_annotated-alles-with-read-support.csv"
    safe_mv_link "$src" "$dst"
done

# Move digger read support results (if they exist)
if [ -d "${orig_outdir}/digger_read_support" ]; then
    for rs_file in ${orig_outdir}/digger_read_support/*_digger_read_support.csv; do
        [ -f "$rs_file" ] || continue
        locus_name=$(basename "$rs_file" _digger_read_support.csv)
        mv "$rs_file" "${outdir}/alleles/${sample}_${locus_name}_digger_annotated-alleles-with-read-support.csv"
        ln -s "${outdir}/alleles/${sample}_${locus_name}_digger_annotated-alleles-with-read-support.csv" "$rs_file"
    done
fi

# Move combined allele tables (if combined mode)
if [[ "$mode" == "combined" ]] && [ -d "${orig_outdir}/combined_alleles" ]; then
    for combined_file in ${orig_outdir}/combined_alleles/*_combined_alleles.csv; do
        [ -f "$combined_file" ] || continue
        fname=$(basename "$combined_file")
        mv "$combined_file" "${outdir}/alleles/${fname}"
        ln -s "${outdir}/alleles/${fname}" "$combined_file"
    done
fi



# Move and link for the contigs-bcftools files
safe_mv_link "${orig_outdir}/vcfs/${sample}_contigs-bcftools_annotated.vcf.gz" "${outdir}/variants/${sample}_contigs-bcftools_annotated.vcf.gz"
safe_mv_link "${orig_outdir}/vcfs/${sample}_contigs-bcftools_annotated.vcf.gz.csi" "${outdir}/variants/${sample}_contigs-bcftools_annotated.vcf.gz.csi"

# Move and link for the ccs-bcftools files
safe_mv_link "${orig_outdir}/vcfs/${sample}_ccs-bcftools_annotated.vcf.gz" "${outdir}/variants/${sample}_ccs-bcftools_annotated.vcf.gz"
safe_mv_link "${orig_outdir}/vcfs/${sample}_ccs-bcftools_annotated.vcf.gz.csi" "${outdir}/variants/${sample}_ccs-bcftools_annotated.vcf.gz.csi"

# Move and link for the clair-from-contigs files
# mv ${orig_outdir}/vcfs/${sample}_clair-from-contigs_annotated.vcf.gz ${outdir}/variants/${sample}_clair-from-contigs_annotated.vcf.gz
# ln -s ${outdir}/variants/${sample}_clair-from-contigs_annotated.vcf.gz ${orig_outdir}/vcfs/${sample}_clair-from-contigs_annotated.vcf.gz

# mv ${orig_outdir}/vcfs/${sample}_clair-from-contigs_annotated.vcf.gz.csi ${outdir}/variants/${sample}_clair-from-contigs_annotated.vcf.gz.csi
# ln -s ${outdir}/variants/${sample}_clair-from-contigs_annotated.vcf.gz.csi ${orig_outdir}/vcfs/${sample}_clair-from-contigs_annotated.vcf.gz.csi

# # Move and link for the clair-from-ccs files
# mv ${orig_outdir}/vcfs/${sample}_clair-from-ccs_annotated.vcf.gz ${outdir}/variants/${sample}_clair-from-ccs_annotated.vcf.gz
# ln -s ${outdir}/variants/${sample}_clair-from-ccs_annotated.vcf.gz ${orig_outdir}/vcfs/${sample}_clair-from-ccs_annotated.vcf.gz

# mv ${orig_outdir}/vcfs/${sample}_clair-from-ccs_annotated.vcf.gz.csi ${outdir}/variants/${sample}_clair-from-ccs_annotated.vcf.gz.csi
# ln -s ${outdir}/variants/${sample}_clair-from-ccs_annotated.vcf.gz.csi ${orig_outdir}/vcfs/${sample}_clair-from-ccs_annotated.vcf.gz.csi




#mv ${orig_outdir}/ccs_cov/average_chrom_coverage.tsv ${outdir}/stats/${sample}_personal-ref-based_depth.tsv
#ln -s ${outdir}/stats/${sample}_personal-ref-based_depth.tsv ${orig_outdir}/ccs_cov/average_chrom_coverage.tsv

safe_mv_link "${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm.stats" "${outdir}/stats/${sample}.asm.stats"
safe_mv_link "${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats" "${outdir}/stats/${sample}.asm-to-ref.flagstats"

safe_mv_link "${orig_outdir}/ccs_cov/${sample}.per-base.bed.gz" "${outdir}/stats/${sample}_ccs_to_ref-based_per-base-depth.bed.gz"
safe_mv_link "${orig_outdir}/ccs_cov/${sample}.per-base.bed.gz.csi" "${outdir}/stats/${sample}_ccs_to_ref-based_per-base-depth.bed.gz.csi"

safe_mv_link "${orig_outdir}/ccs_cov/${sample}.regions.bed.gz" "${outdir}/stats/${sample}_ccs_to_ref-based_regions-depth.bed.gz"
safe_mv_link "${orig_outdir}/ccs_cov/${sample}.regions.bed.gz.csi" "${outdir}/stats/${sample}_ccs_to_ref-based_regions-depth.bed.gz.csi"

safe_mv_link "${orig_outdir}/${sample}_readLengthHistogram.png" "${outdir}/stats/${sample}_readLengthHistogram.png"

pigz -p "${threads}" ${outdir}/reads/ccs-reads.fasta