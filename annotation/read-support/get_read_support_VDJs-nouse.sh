#!/bin/bash
set -x

reffn=$2
IG_loci=$3
threads=$4
#masked_ref=/home/egenge01/projects/12_sample_test/reference_IGloci_masked.fasta
scratch=$5
minimap_option=$6
mask_ref=${scratch}/ref_IG_masked.fasta

function run_make_ref_masked {

bedtools maskfasta -fi ${reffn} -bed ${IG_loci} -fo ${scratch}/ref_IG_masked.fasta
samtools faidx ${scratch}/ref_IG_masked.fasta

}

function run_map_ccs_to_pers {
    while IFS=$'\t' read -r sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import trb_gene trb_import trg_gene trg_import trd_gene trd_import tra_gene tra_import ccs_bam
    do
    echo sample is ${sample}
    echo ccs bam is ${ccs_bam}
	mkdir -p ${scratch}/read_support/${sample}
	outd=${scratch}/read_support/${sample}
	mkdir -p ${outd}/ccs_to_pers
	#convert pacbio hifi reads to fasta
	#if [[ "${ccs_bam}" == *.bam ]]; then
    #    samtools view ${ccs_bam} | awk '{ print ">"$1"\n"$10 }' > ${outd}/ccs_to_pers/reads.fasta
    #elif [[ "${ccs_bam}" == *.fasta ]]; then
    #    cp ${ccs_bam} ${outd}/ccs_to_pers/reads.fasta
    #else
    #    echo "Unsupported file format"
    #fi

	samtools faidx ${scratch}/reads.fasta
	
	#create personalized reference
	samtools view -F 0x100 -F 0x800 "${asm_bam}" | awk '{print ">"$1"\n"$10}' > "${outd}/ccs_to_pers/contigs.fasta"
	samtools faidx ${outd}/ccs_to_pers/contigs.fasta
	pers_contigs=${outd}/ccs_to_pers/contigs.fasta
	cat ${mask_ref} ${pers_contigs} \
	    > ${outd}/ccs_to_pers/pers_ref.fasta
	samtools faidx ${outd}/ccs_to_pers/pers_ref.fasta
	#make not gpu one day
	minimap2 -a ${minimap_option} --secondary=yes -t ${threads} -L ${outd}/ccs_to_pers/pers_ref.fasta ${scratch}/reads.fasta > ${outd}/ccs_to_pers/output.sam;
    samtools view -Sbh ${outd}/ccs_to_pers/output.sam > ${outd}/ccs_to_pers/output.bam;
    samtools sort -@ ${threads} ${outd}/ccs_to_pers/output.bam -o ${outd}/ccs_to_pers/output.sorted.bam;
    samtools index ${outd}/ccs_to_pers/output.sorted.bam;
    rm -f ${outd}/ccs_to_pers/output.sam
    done < $fofn
}

function run_append_pos {
    while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import trb_gene trb_import trg_gene trg_import trd_gene trd_import tra_gene tra_import ccs_bam
    do
        # Base directory for imported genes
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        
        # Create subdirectories for each type and define output paths
        mkdir -p "${base_outd}/IGK"
        chr2_import_out="${base_outd}/IGK/$(basename "${chr2_import}")"

        mkdir -p "${base_outd}/IGL"
        chr22_import_out="${base_outd}/IGL/$(basename "${chr22_import}")"

        mkdir -p "${base_outd}/IGH"
        igh_import_out="${base_outd}/IGH/$(basename "${igh_import}")"

        mkdir -p "${base_outd}/IGHC"
        ighc_import_out="${base_outd}/IGHC/$(basename "${ighc_import}")"

        mkdir -p "${base_outd}/TRB"
        trb_import_out="${base_outd}/TRB/$(basename "${trb_import}")"

        mkdir -p "${base_outd}/TRG"
        trg_import_out="${base_outd}/TRG/$(basename "${trg_import}")"

        mkdir -p "${base_outd}/TRD"
        trd_import_out="${base_outd}/TRD/$(basename "${trd_import}")"
       
        mkdir -p "${base_outd}/TRA"
        tra_import_out="${base_outd}/TRA/$(basename "${tra_import}")"

        #mkdir -p "${base_outd}/TRA"
        #ighc_import_out="${base_outd}/TRA/$(basename "${tra_import}")"

        # Call the /opt/wasp/conda/bin/python script with input, original output, and modified output paths
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${chr2_gene}" "${chr2_import}" "${chr2_import_out}"
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${chr22_gene}" "${chr22_import}" "${chr22_import_out}"
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${igh_gene}" "${igh_import}" "${igh_import_out}"
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${trb_gene}" "${trb_import}" "${trb_import_out}"
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${trg_gene}" "${trg_import}" "${trg_import_out}"
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${tra_gene}" "${tra_import}" "${tra_import_out}"
        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/append_pos_import_genes.py "${trd_gene}" "${trd_import}" "${trd_import_out}"

        /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/ighc_append_pos.py "${ighc_gene}" "${ighc_import}" "${ighc_import_out}"
#        /opt/wasp/conda/bin/python append_pos_import_genes_ighc.py "${ighc_gene}" "${ighc_import}" "${ighc_import_out}"
    done < $fofn
}

function get_read_support_vdj3 {
    while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import trb_gene trb_import trg_gene trg_import trd_gene trd_import tra_gene tra_import ccs_bam
    do
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        bam_file="${scratch}/read_support/${sample}/ccs_to_pers/output.sorted.bam" # ccs reads to personalized reference
        ref="${scratch}/read_support/${sample}/ccs_to_pers/pers_ref.fasta" # personalized reference

        if [ ! -f "${bam_file}.bai" ]; then
            samtools index "$bam_file"
        fi

        #for gene_type in "chr2" "chr22" "igh" "trb" "trg" "tra" "trd" #"ighc"
        for gene_type in "IGK" "IGL" "IGH" "TRB" "TRG" "TRD" "TRA" #"ighc"
        do
            import_out="${base_outd}/${gene_type}/${sample}_make_gene_file_imported.csv"

            if [[ -f "$import_out" ]]; then
		modified_import_out="$import_out"
		
                tmp_file="${import_out}_read_support.tmp"
                echo "Total_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Percent_Accuracy,Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match" > "$tmp_file"
                header=$(head -n 1 "$import_out")
                IFS=',' read -ra header_cols <<< "$header"
                for i in "${!header_cols[@]}"; do
                    case "${header_cols[$i]}" in
                        "contig") contig_col=$i ;;
                        "REGION_start") start_col=$i ;;
                        "REGION_end") end_col=$i ;;
			"gene") gene_col=$i ;;  # Added case for 'gene'
                    esac
                done
		tmp_counts="${tmp_file}_counts"  # Temporary file to hold counts
		> "$tmp_counts"  # Clear or create the temp file for counts

                tail -n +2 "$import_out" | while IFS=, read -ra line
                do
                    contig="${line[$contig_col]}"
                    start=$(echo "${line[$start_col]}" | awk '{printf "%.0f", $1}')
                    end=$(echo "${line[$end_col]}" | awk '{printf "%.0f", $1}')
                    
		    gene="${line[$gene_col]}"  # Extract the 'gene' value using the identified column index
                    region="${contig}:${start}-${end}"

                    contig_filename=$(echo "$contig" | tr '/' '_')
                    tmp_bam="${base_outd}/${gene_type}/${contig_filename}_${start}_${end}.bam"
                    mkdir -p "$(dirname "$tmp_bam")"
                    echo start is $start. end is $end. contig is $contig. 
                    samtools view -F 0x100 -F 0x800 -b "$bam_file" -o "$tmp_bam" -U "/dev/null" "${contig}:${start}-${end}"
                    samtools index "$tmp_bam"

		    samtools mpileup -f "$ref" -r "$region" "$tmp_bam" | \
			awk -v total_positions="$((end - start + 1))" -v sample="$sample" \
			'BEGIN {
    total_reads=0; mismatched_positions=0; matched_positions=0; positions_with_10x=0;
    mismatch_list=""; match_list="";
}
{
    total_reads += length($5);
    mismatches = length(gensub(/[.,]/, "", "g", $5));
    matches = length(gensub(/[^.,]/, "", "g", $5));
    mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches);
    match_list = (match_list == "" ? matches : match_list ":" matches);

    coverage = length($5);
    if (coverage >= 10) {
        positions_with_10x++;
    }

    mismatch_rate = mismatches / coverage;
    match_rate = matches / coverage;

    if (mismatch_rate > 0.2) {
        mismatched_positions++;
    }

    if (match_rate > 0.8) {
        matched_positions++;
    }
}
END {
    avg_reads_per_position = (total_positions > 0) ? total_reads / total_positions : 0;
    percent_accuracy = (matched_positions / total_positions) * 100;
    print total_positions, avg_reads_per_position, mismatched_positions, matched_positions, mismatch_list, match_list, percent_accuracy, positions_with_10x;}' OFS=',' >> "${tmp_file}_awk_out"
		    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/match_subsequences.py "$tmp_bam" "$contig" "$start" "$end" "$gene" "$import_out" > "${tmp_file}_py_out"
		    wait
		    paste -d ',' "${tmp_file}_awk_out" "${tmp_file}_py_out" >> "$tmp_file"
		    rm "${tmp_file}_awk_out" "${tmp_file}_py_out"
                    rm "$tmp_bam" "${tmp_bam}.bai"
                done

# New block to merge data and update Subject and Sample_Name
		if [[ -f "$import_out" && -f "$tmp_file" ]]; then
		    combined_file="${import_out%.csv}_combined.csv"
		    final_output="${import_out%.csv}_with_read_support.csv"
		    
		    # Merge the original import file with the new tmp file
		    paste -d ',' "$import_out" "$tmp_file" > "$combined_file"
		    rm -f "$tmp_file"
		    
		    # Update the Subject and Sample_Name columns in the final output
		    awk -v sample="$sample" 'BEGIN{FS=OFS=","} {
    if (NR == 1) {
        print;
    } else {
        $2 = sample;  # Update the 2nd column with $sample
        $3 = sample;  # Update the 3rd column with $sample
        print;
    }
}' "$combined_file" > "$final_output"

		    rm -f "$combined_file"
		fi
	    fi
        done
    done < $fofn
}


function get_read_support_ighc {
cat $fofn | while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import trb_gene trb_import trg_gene trg_import trd_gene trd_import tra_gene tra_import ccs_bam
    do
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        bam_file="${scratch}/read_support/${sample}/ccs_to_pers/output.sorted.bam" # ccs reads to personalized reference
        ref="${scratch}/read_support/${sample}/ccs_to_pers/pers_ref.fasta" # personalized reference

        if [ ! -f "${bam_file}.bai" ]; then
            samtools index "$bam_file"
        fi

        for gene_type in "IGHC" #"chr2" "chr22" "igh"
        do
             import_out="${base_outd}/${gene_type}/${sample}_make_gene_file_imported.csv"
             #import_out="${base_outd}/${gene_type}/IGenotyper_imported_genes.csv"
            if [[ -f "$import_out" ]]; then
                # Use csvcut to remove the "notes" column, which causes problems and can have or always has a comma
                #csvcut -C "notes" "$import_out" > "${import_out%.csv}_nonotes.csv"
                #modified_import_out="${import_out%.csv}_nonotes.csv"


                modified_import_out=$import_out


                tmp_file="${import_out}_read_support.tmp"

		echo "Total_Positions,Total_Reads_by_Positions,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Mismatched_Positions_Coverage_Less_Than_10,Mismatched_Positions_Coverage_10_Or_Greater,Matched_Positions_Coverage_Less_Than_10,Matched_Positions_Coverage_10_Or_Greater,Percent_Accuracy,Fully_Spanning_Allele_reads,Fully_Spanning_Allele_reads_100_Match,Allele_reads_100_Match_e1,Allele_reads_100_Match_e2,Allele_reads_100_Match_e3,Allele_reads_100_Match_e4,Allele_reads_100_Match_e5,Allele_reads_100_Match_e6,Allele_reads_100_Match_e7,Allele_reads_100_Match_e8,Allele_reads_100_Match_e9" > "$tmp_file"

		header=$(head -n 1 "$modified_import_out")
		IFS=',' read -ra header_cols <<< "$header"
		declare -A exon_cols  # Use an associative array to store exon column indices
		
		for i in "${!header_cols[@]}"; do
		    case "${header_cols[$i]}" in
			"contig") contig_col=$i ;;
			C-EXON_[1-9]_start) exon_cols["${header_cols[$i]}"]=$i ;;  # Match any C-EXON_X_start where X is 1-9
			C-EXON_[1-9]_end) exon_cols["${header_cols[$i]}"]=$i ;;  # Match any C-EXON_X_end where X is 1-9
			"gene") gene_col=$i ;;  # Added case for 'gene'
		    esac
		done
		
		# To access a specific exon start or end column index later in the script, you can use:
		# ${exon_cols[C-EXON_4_start]} or ${exon_cols[C-EXON_4_end]}, etc.
		
		# Specify the output file name
#		output_file="exon_regions.tsv"
		
		# Write the header to the output file
	#	echo -e "Contig\tExon_Start\tExon_End" 
#		> "$output_file"
		
		# Iterate over each line of the input file, skipping the header
		tail -n +2 "$modified_import_out" | while IFS=, read -ra line; do
		    contig="${line[$contig_col]}"
		    
		    # Generate a unique filename for the BED file based on the contig name
		    contig_filename=$(echo "$contig" | tr '/' '_')
		    bed_file="${base_outd}/${gene_type}/${contig_filename}.bed"
		    mkdir -p "$(dirname "$bed_file")"
		    
		    # Clear the BED file for new entries
		    > "$bed_file"
		    
		    # Iterate over each exon start and end index pair from the exon_cols associative array
# Iterate over each exon start and end index pair from the exon_cols associative array
            for key in "${!exon_cols[@]}"; do
                if [[ $key =~ C-EXON_[1-9]_start ]]; then
                    exon_start=$(echo "${line[${exon_cols[$key]}]}" | awk '{printf "%.0f", $1-1}')
                    exon_end_key="${key/_start/_end}"  # Replace 'start' with 'end' to get the corresponding end key
                    exon_end=$(echo "${line[${exon_cols[$exon_end_key]}]}" | awk '{printf "%.0f", $1}')

                            # Write the contig, exon start, and exon end to the BED file
                    echo -e "${contig}\t${exon_start}\t${exon_end}" | awk '{if ($3 !=0) print $0}' >> "$bed_file"
                fi
            done
		    
		    gene="${line[$gene_col]}"  # Extract the 'gene' value using the identified column index
		    

#                    contig_filename=$(echo "$contig" | tr '/' '_')
                    tmp_bam="${base_outd}/${gene_type}/${contig_filename}_${start}_${end}.bam"
                    mkdir -p "$(dirname "$tmp_bam")"
                    echo start is $start. end is $end. contig is $contig. 
                    samtools view -F 0x100 -F 0x800 -b "$bam_file" -o "$tmp_bam" -U "/dev/null" "${contig}:${start}-${end}"
                    samtools index "$tmp_bam"
                    
		    total_positions=$(awk '{sum += ($3 - $2)} END {print sum}' "$bed_file")
                    samtools mpileup -f "$ref" -l "$bed_file" "$tmp_bam" | \
			awk -v total_positions="$total_positions" \
                    'BEGIN {
                        total_reads=0; mismatched_positions=0; matched_positions=0;
                        mismatch_list=""; match_list="";
                        mismatched_positions_coverage_less_than_10=0;
                        mismatched_positions_coverage_10_or_greater=0;
                        matched_positions_coverage_less_than_10=0;
                        matched_positions_coverage_10_or_greater=0;
                    }
                    {
                        total_reads += $4;
                        #mismatches = length(gensub(/[.,]/, "", "g", $5));
                        matches = length(gensub(/[^.,]/, "", "g", $5));
                        mismatches = $4 - matches
                        mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches);
                        match_list = (match_list == "" ? matches : match_list ":" matches);

                        coverage = $4;
                        mismatch_rate = mismatches / coverage;
                        match_rate = matches / coverage;

                        if (mismatch_rate > 0.2) {
                            mismatched_positions++;
                            if (coverage < 10) {
                                mismatched_positions_coverage_less_than_10++;
                            } else {
                                mismatched_positions_coverage_10_or_greater++;
                            }
                        }

                        if (match_rate >= 0.8) {
                            matched_positions++;
                            if (coverage < 10) {
                                matched_positions_coverage_less_than_10++;
                            } else {
                                matched_positions_coverage_10_or_greater++;
                            }
                        }
                    }
                    END {
                        percent_accuracy = (matched_positions / total_positions) * 100;
                        print total_positions, total_reads, mismatched_positions, matched_positions, mismatch_list, match_list, mismatched_positions_coverage_less_than_10, mismatched_positions_coverage_10_or_greater, matched_positions_coverage_less_than_10, matched_positions_coverage_10_or_greater, percent_accuracy;
                    }' OFS=',' >> "${tmp_file}_awk_out"
		    /opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/read-support/ighc_match3.py "$tmp_bam" "$contig" "$gene" "$modified_import_out" > "${tmp_file}_py_out"
		    wait
		    paste -d ',' "${tmp_file}_awk_out" "${tmp_file}_py_out" >> "$tmp_file"
		    rm "${tmp_file}_awk_out" "${tmp_file}_py_out"
		    rm "$tmp_bam" "${tmp_bam}.bai"
		done
#		New block to merge data and update Subject and Sample_Name
		if [[ -f "$modified_import_out" && -f "$tmp_file" ]]; then
		    combined_file="${modified_import_out%.csv}_combined.csv"
		    final_output="${modified_import_out%.csv}_with_read_support.csv"
		    
		    # Merge the original import file with the new tmp file
		    paste -d ',' "$modified_import_out" "$tmp_file" > "$combined_file"
		    rm -f "$tmp_file"
		    
		    # Update the Subject and Sample_Name columns in the final output
		    awk -v sample="$sample" 'BEGIN{FS=OFS=","} {
    if (NR == 1) {
        print;
    } else {
        $2 = sample;  # Update the 2nd column with $sample
        $3 = sample;  # Update the 3rd column with $sample
        print;
    }
}' "$combined_file" > "$final_output"
		    
		    rm -f "$combined_file"
		fi
	    fi
	done
done
}

fofn=$1
run_make_ref_masked
run_map_ccs_to_pers
run_append_pos
get_read_support_vdj3
get_read_support_ighc
rm ${scratch}/ref_IG_masked.fasta
