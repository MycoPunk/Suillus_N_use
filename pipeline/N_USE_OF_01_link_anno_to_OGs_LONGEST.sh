#!/bin/bash
#SBATCH --job-name=link_annos
#SBATCH --output=logs/link_anno_to_OGs.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=1:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=common

#input file
input_file="N_USE_longest_headers_clean.txt"

#output file
output_file="merged_annotations.txt"

#clear the output file if it exists
> "$output_file"

#read the input file line by line
while IFS=' ' read -r og_number species protein_id; do
    #construct the path to the annotation file
    anno_file="genomes/${species}/anno/${species}_GeneCatalog_proteins_*_IPR.tab"
    
    #find all matching lines in the annotation file and process them
    grep "^${protein_id}" $anno_file | while IFS=$'\t' read -r _ annotation_rest; do
        echo -e "${og_number}\t${species}\t${protein_id}\t${annotation_rest}" >> "$output_file"
    done
    
    #if no matching lines found, print an error message
    if [ $? -ne 0 ]; then
        echo "No annotation found for ${og_number} ${species} ${protein_id}" >&2
    fi
done < "$input_file"

echo "Processing complete. Results written to ${output_file}"
