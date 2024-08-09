#!/bin/bash
#SBATCH --job-name=make_merops_anno_map
#SBATCH --output=logs/make_merops_anno_map.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=1:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=common

# Input file path
input_file="/hpc/group/vilgalyslab/lal76/databases/merops_scan_peptidase_only.lib"

# Output file name
output_file="merops_anno_map.txt"

# Process the file
awk '
/^>/ {
    #extract the first word after ">"
    first_word = substr($1, 2)  # Remove the ">" character

    #find the position of " - " and " ("
    start = index($0, " - ") + 3
    end = index($0, " (")

    #extract the text between " - " and " ("
    description = substr($0, start, end - start)

    #print the result
    print first_word "\t" description
}
' "$input_file" > "$output_file"

echo "Processing complete. Results written to $output_file"
