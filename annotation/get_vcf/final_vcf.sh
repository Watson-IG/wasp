#!/bin/bash
set -e

# Define the main directory for VCF files
scratch="$PWD"
outd="$scratch/run_wasp"

function genotype_SV_regions {
    local bam_path="$1"
    local SV_regions_1bp="$2"
    local sample="$3"

    # Ensure the directory exists
    local sample_vcf_dir="${outd}/${sample}/vcfs"
    mkdir -p "$sample_vcf_dir"

    local sample_sv_results="${sample_vcf_dir}/${sample}_sv_genotype_results.txt"

    while read -r sv_region; do
        sv_name=$(echo "$sv_region" | cut -f4)
        grep -w "$sv_name" "$SV_regions_1bp" > "${outd}/${sample}/vcfs/${sv_name}.bed"

        mpileup_output=$(samtools mpileup -l "${outd}/${sample}/vcfs/${sv_name}.bed" -f "${reffn}" "$bam_path" 2>/dev/null | head -1)
        genotype=$(echo "$mpileup_output" | awk '{print $5}')
        asterisk_count=$(echo "$genotype" | tr -cd '*' | wc -c)
        total_count=$(echo "$genotype" | wc -c)

        if [ "$asterisk_count" -eq 0 ]; then
            genotype_label="0/0"
        elif [ "$asterisk_count" -eq "$((total_count-1))" ]; then
            genotype_label="1/1"
        else
            genotype_label="0/1"
        fi

        echo -e "$sample\t$sv_name\t$genotype_label" >> "$sample_sv_results"
    done < "$SV_regions_1bp"
}

function process_vcf {
    local bam_file="$1"
    local sample="$2"
    local reffn="$3"
    local num_threads="$4"
    local SV_regions_entire="$5"
    local changeg="$6"
    local anno_config_file="$7"
    local vcfanno="$8"
    local bed_dir="$9"
    local ccs_bam_file="${10}"
    
    local sample_vcf_dir="${outd}/${sample}/vcfs"
    mkdir -p "$sample_vcf_dir"
    local anno_log="${sample_vcf_dir}/vcf_annotation.log"

    local of="${sample_vcf_dir}/${sample}"
    echo "Running bcftools mpileup and call on contigs..."
    bcftools mpileup -B -a QS -Ou -f "${reffn}" -R ${bed_dir}/IG_loci.bed --threads "${num_threads}" "$bam_file" 2>>"${anno_log}" | \
    bcftools call -m -Oz -o "${of}.tmp.contigs.vcf.gz" 2>>"${anno_log}"
    bcftools sort --write-index "${of}.tmp.contigs.vcf.gz" -Oz -o "${of}.contigs.vcf.gz" 2>>"${anno_log}"

    echo "Running bcftools mpileup and call on CCS reads..."
    bcftools mpileup -B -a QS -Ou -f "${reffn}" -R ${bed_dir}/IG_loci.bed --threads "${num_threads}" "$ccs_bam_file" 2>>"${anno_log}" | \
    bcftools call -m -Oz -o "${of}.tmp.ccs.vcf.gz" 2>>"${anno_log}"
    bcftools sort --write-index "${of}.tmp.ccs.vcf.gz" -Oz -o "${of}.ccs.vcf.gz" 2>>"${anno_log}"

    echo "Annotating contigs VCF with vcfanno..."
    "${vcfanno}" "${anno_config_file}" "${of}.contigs.vcf.gz" 2>>"${anno_log}" > "${sample_vcf_dir}/${sample}_contigs-bcftools_annotated.vcf"
    bgzip -f "${sample_vcf_dir}/${sample}_contigs-bcftools_annotated.vcf"
    bcftools index "${sample_vcf_dir}/${sample}_contigs-bcftools_annotated.vcf.gz" 2>>"${anno_log}"

    echo "Annotating CCS VCF with vcfanno..."
    "${vcfanno}" "${anno_config_file}" "${of}.ccs.vcf.gz" 2>>"${anno_log}" > "${sample_vcf_dir}/${sample}_ccs-bcftools_annotated.vcf"
    bgzip -f "${sample_vcf_dir}/${sample}_ccs-bcftools_annotated.vcf"
    bcftools index "${sample_vcf_dir}/${sample}_ccs-bcftools_annotated.vcf.gz" 2>>"${anno_log}"
}

# Inputs from the user or script
sample="$1"
bam_file="$2"
reffn="$3"
num_threads="$4"
SV_regions_entire="/opt/wasp/scripts/annotation/KL_SV_regions_entire.bed"
SV_regions_1bp="/opt/wasp/scripts/annotation/SV_regions_1bp.bed"
changeg="/opt/wasp/scripts/annotation/get_vcf/vcf_processing.py"
anno_config_file="/opt/wasp/scripts/annotation/config.toml"
vcfanno="vcfanno"
bed_dir=$5
ccs_bam_file="$6"

# Add and index new read group
samtools addreplacerg -r ID:"${sample}" -r SM:"${sample}" -o "${outd}/$sample/${sample}.editRG.bam" "${bam_file}" 2>/dev/null
samtools index "${outd}/$sample/${sample}.editRG.bam" 2>/dev/null

# Run functions
genotype_SV_regions "${outd}/$sample/${sample}.editRG.bam" "$SV_regions_1bp" "$sample"
process_vcf "${outd}/$sample/${sample}.editRG.bam" "$sample" "$reffn" "$num_threads" "$SV_regions_entire" "$changeg" "$anno_config_file" "$vcfanno" "$bed_dir" "$ccs_bam_file"
