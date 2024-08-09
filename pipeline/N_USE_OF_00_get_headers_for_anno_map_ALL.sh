#!/bin/bash
#SBATCH --job-name=get_headers_for_anno_map   
#SBATCH --output=logs/get_headers_for_anno_map.%A.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --partition=common

#create the output directory if it doesn't exist
mkdir -p headers

#process each .fa file in the specified directory
for input_file in /hpc/group/vilgalyslab/lal76/Suillus/N_use/OF_results/Results_May22/Orthogroup_Sequences/*.fa
do
    #extract the base name of the input file (this will be the OG number)
    base_name=$(basename "$input_file" .fa)
    
    #define the output file name
    output_file="headers/${base_name}_headers.txt"
    
    #process the file
    awk -v og="$base_name" '
    /^>/ {
        sub(/^>/, "");  # Remove the ">" at the start
        split($0, parts, "|");  # Split the header on "|"
        split(parts[2], id_parts, "_");  # Split the second part on "_"
        printf "%s\t%s\t%s\n", og, parts[1], id_parts[2];
    }
    ' "$input_file" > "$output_file"
    
    echo "Processed $input_file -> $output_file"
done

echo "All files processed. Output files are in the 'headers' directory."
