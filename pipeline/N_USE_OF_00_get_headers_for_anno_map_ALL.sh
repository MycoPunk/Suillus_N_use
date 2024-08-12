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

##this script pulls out the headers from the fastas of each OG, creating one output file for each orthogroup with 4 cols
##1) the OG#, the sp, the protin ID, and the length of each protein sequence 

#create the output directory if it doesn't exist
mkdir -p headers

#process each .fa file in the OG directory
for input_file in /hpc/group/vilgalyslab/lal76/Suillus/N_use/OF_results/Results_May22/Orthogroup_Sequences/*.fa
do
    #extract the base name of the input file (this will be the OG number)
    base_name=$(basename "$input_file" .fa)
    
    #define the output file name
    output_file="headers/${base_name}_headers.txt"
    
    #process the file
    sed -e 's/\r//g' "$input_file" | awk -v og="$base_name" '
    BEGIN {
        RS = ">"
        FS = "\n"
    }
    NR > 1 {
        header = $1
        split(header, parts, "|")
        species = parts[1]
        split(parts[2], id_parts, "_")
        id = id_parts[2]
        
        seq = ""
        for (i=2; i<=NF; i++) {
            seq = seq $i
        }
        gsub(/[^A-Za-z]/, "", seq)
        
        printf "%s\t%s\t%s\t%d\n", og, species, id, length(seq)
    }
    ' > "$output_file"
    
    echo "Processed $input_file -> $output_file"
done

echo "All files processed. Output files are in the 'headers' directory."
