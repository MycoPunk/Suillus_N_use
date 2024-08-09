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

#this script links the interpro annotations to the OGs, creating a seprate output file for each orthogroup annotations for each representative.
#where no annotations are available, it creates an empty file. 

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
    
    #read the input file line by line
    while IFS=$'\t' read -r og_number species protein_id; do
        # Construct the path to the annotation file
        anno_file="genomes/${species}/anno/${species}_GeneCatalog_proteins_*_IPR.tab"
        
        #find all matching lines in the annotation file and process them
        grep "^${protein_id}" $anno_file | while IFS=$'\t' read -r _ annotation_rest; do
            echo -e "${og_number}\t${species}\t${protein_id}\t${annotation_rest}" >> "$output_file"
        done
        
        #if no matching lines are found, print an error message
        if [ $? -ne 0 ]; then
            echo "No annotation found for ${og_number} ${species} ${protein_id}" >&2
        fi
    done < "$input_file"
    
    echo "Processing complete for ${input_file}. Results written to ${output_file}"
done

echo "All files processed. Output files are in the 'anno' directory."
