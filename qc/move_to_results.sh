#!/bin/bash

sample=$1
orig_outdir=$2
outdir=$PWD/results/${sample}

mkdir -p ${outdir}
mkdir -p ${outdir}/reads ${outdir}/alignments ${outdir}/variants ${outdir}/alleles ${outdir}/stats

# Moving and creating symlinks
mv ${orig_outdir}/reads.fasta ${outdir}/reads/ccs-reads.fasta
ln -s ${outdir}/reads/ccs-reads.fasta ${orig_outdir}/reads.fasta


mv ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta ${outdir}/reads/hifiasm_ig-filtered_contigs.fasta
ln -s ${outdir}/reads/hifiasm_ig-filtered_contigs.fasta ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/contigs.fasta

mv ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${outdir}/alignments/${sample}_contigs-to-ref.sorted.bam
ln -s ${outdir}/alignments/${sample}_contigs-to-ref.sorted.bam ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam

mv ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam.bai ${outdir}/alignments/${sample}_contigs-to-ref.sorted.bam.bai
ln -s ${outdir}/alignments/${sample}_contigs-to-ref.sorted.bam.bai ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam.bai

mv ${orig_outdir}/ccs_cov/ccs_to_ref.sorted.bam ${outdir}/alignments/${sample}_ccs-to-ref.sorted.bam
ln -s ${outdir}/alignments/${sample}_ccs-to-ref.sorted.bam ${orig_outdir}/ccs_cov/ccs_to_ref.sorted.bam

mv ${orig_outdir}/ccs_cov/ccs_to_ref.sorted.bam.bai ${outdir}/alignments/${sample}_ccs-to-ref.sorted.bam.bai
ln -s ${outdir}/alignments/${sample}_ccs-to-ref.sorted.bam.bai ${orig_outdir}/ccs_cov/ccs_to_ref.sorted.bam.bai

mv ${orig_outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam ${outdir}/alignments/${sample}_ccs-to-personal-reference.sorted.bam
ln -s ${outdir}/alignments/${sample}_ccs-to-personal-reference.sorted.bam ${orig_outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam

#mv ${orig_outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam.bai ${outdir}/alignments/${sample}_ccs-to-personal-reference.sorted.bam.bai
#ln -s ${outdir}/alignments/${sample}_ccs-to-personal-reference.sorted.bam.bai ${orig_outdir}/read_support/${sample}/ccs_to_pers/output.sorted.bam.bai

# Loci list
loci_list=("IGH" "IGHC" "IGK" "IGL" "TRA" "TRB" "TRD" "TRG")

for loci in "${loci_list[@]}"; do
    #mkdir -p "${outdir}/alleles/${loci}"
    mv "${orig_outdir}/read_support/${sample}/imported_genes/${loci}/${sample}_make_gene_file_imported_with_read_support.csv" \
       "${outdir}/alleles/${sample}_${loci}_annotated-alles-with-read-support.csv"
    ln -s "${outdir}/alleles/${sample}_${loci}_annotated-alles-with-read-support.csv" \
          "${orig_outdir}/read_support/${sample}/imported_genes/${loci}/${sample}_make_gene_file_imported_with_read_support.csv"
done



# Move and link for the contigs-bcftools files
mv ${orig_outdir}/vcfs/${sample}_contigs-bcftools_annotated.vcf.gz ${outdir}/variants/${sample}_contigs-bcftools_annotated.vcf.gz
ln -s ${outdir}/variants/${sample}_contigs-bcftools_annotated.vcf.gz ${orig_outdir}/vcfs/${sample}_contigs-bcftools_annotated.vcf.gz

mv ${orig_outdir}/vcfs/${sample}_contigs-bcftools_annotated.vcf.gz.csi ${outdir}/variants/${sample}_contigs-bcftools_annotated.vcf.gz.csi
ln -s ${outdir}/variants/${sample}_contigs-bcftools_annotated.vcf.gz.csi ${orig_outdir}/vcfs/${sample}_contigs-bcftools_annotated.vcf.gz.csi

# Move and link for the ccs-bcftools files
mv ${orig_outdir}/vcfs/${sample}_ccs-bcftools_annotated.vcf.gz ${outdir}/variants/${sample}_ccs-bcftools_annotated.vcf.gz
ln -s ${outdir}/variants/${sample}_ccs-bcftools_annotated.vcf.gz ${orig_outdir}/vcfs/${sample}_ccs-bcftools_annotated.vcf.gz

mv ${orig_outdir}/vcfs/${sample}_ccs-bcftools_annotated.vcf.gz.csi ${outdir}/variants/${sample}_ccs-bcftools_annotated.vcf.gz.csi
ln -s ${outdir}/variants/${sample}_ccs-bcftools_annotated.vcf.gz.csi ${orig_outdir}/vcfs/${sample}_ccs-bcftools_annotated.vcf.gz.csi

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

mv ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm.stats ${outdir}/stats/${sample}.asm.stats
ln -s ${outdir}/stats/${sample}.asm.stats ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm.stats

mv ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats ${outdir}/stats/${sample}.asm-to-ref.flagstats
ln -s ${outdir}/stats/${sample}.asm-to-ref.flagstats ${orig_outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats

mv ${orig_outdir}/ccs_cov/${sample}.per-base.bed.gz ${outdir}/stats/${sample}_ccs_to_ref-based_per-base-depth.bed.gz
ln -s ${outdir}/stats/${sample}_ccs_to_ref-based_per-base-depth.bed.gz ${orig_outdir}/ccs_cov/${sample}.per-base.bed.gz
mv ${orig_outdir}/ccs_cov/${sample}.per-base.bed.gz.csi ${outdir}/stats/${sample}_ccs_to_ref-based_per-base-depth.bed.gz.csi
ln -s ${outdir}/stats/${sample}_ccs_to_ref-based_per-base-depth.bed.gz.csi ${orig_outdir}/ccs_cov/${sample}.per-base.bed.gz.csi

mv ${orig_outdir}/ccs_cov/${sample}.regions.bed.gz ${outdir}/stats/${sample}_ccs_to_ref-based_regions-depth.bed.gz
ln -s ${outdir}/stats/${sample}_ccs_to_ref-based_regions-depth.bed.gz ${orig_outdir}/ccs_cov/${sample}.regions.bed.gz
mv ${orig_outdir}/ccs_cov/${sample}.regions.bed.gz.csi ${outdir}/stats/${sample}_ccs_to_ref-based_regions-depth.bed.gz.csi
ln -s ${outdir}/stats/${sample}_ccs_to_ref-based_regions-depth.bed.gz.csi ${orig_outdir}/ccs_cov/${sample}.regions.bed.gz.csi

mv ${orig_outdir}/${sample}_readLengthHistogram.png ${outdir}/stats/${sample}_readLengthHistogram.png
ln -s ${outdir}/stats/${sample}_readLengthHistogram.png ${orig_outdir}/${sample}_readLengthHistogram.png

gzip ${outdir}/reads/ccs-reads.fasta