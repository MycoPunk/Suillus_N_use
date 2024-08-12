#!/bin/bash
#SBATCH --job-name=link_annos
#SBATCH --output=logs/link_anno_to_OGs.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=18:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=common

#this script links the interpro annotations to the OGs, creating a separate output file for each orthogroup. It privelages the longest protein file to base the annotation on
#based on the last column of the _headers.txt files generated in the previous step of the pipeline. If no annotation is found for the longest OG representative, 
#the script continues by looking for an annotation for the second longest protein, and so on. 
#where no annotations are available for any members of the OG, it creates an empty file.


#create the output directory if it doesn't exist
mkdir -p anno

#process each _headers.txt file in the headers directory
for input_file in headers/*_headers.txt
do
    #extract the base name of the input file (without _headers.txt)
    base_name=$(basename "$input_file" _headers.txt)
    
    #define the output file name
    output_file="anno/${base_name}_anno.txt"
    
    #clear the output file if it exists
    > "$output_file"
    
    #sort the input file by the fourth column (length) in descending order
    sorted_input=$(sort -k4 -nr "$input_file")
    
    #process each line of the sorted input
    echo "$sorted_input" | while IFS=$'\t' read -r og_number species protein_id length; do
        # Construct the path to the annotation file
        anno_file="genomes/${species}/anno/${species}_GeneCatalog_proteins_*_IPR.tab"
        
        #find all matching lines in the annotation file and process them
        awk -v pid="$protein_id" '$1 == pid' $anno_file | while IFS=$'\t' read -r _ annotation_rest; do
            echo -e "${og_number}\t${species}\t${protein_id}\t${annotation_rest}" >> "$output_file"
        done
        
        #if annotations were found and written to the output file, break the loop
        if [ -s "$output_file" ]; then
            break
        fi
        
        #if no matching lines are found, print an error message
        if [ $? -ne 0 ]; then
            echo "No annotation found for ${og_number} ${species} ${protein_id}" >&2
        fi
    done
    
    #if the output file is empty after processing all proteins, print a warning
    if [ ! -s "$output_file" ]; then
        echo "Warning: No annotations found for any proteins in ${base_name}" >&2
    else
        echo "Processing complete for ${input_file}. Results written to ${output_file}"
    fi
done

echo "All files processed. Output files are in the 'anno' directory."
