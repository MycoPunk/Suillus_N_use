#!/bin/bash
#SBATCH --job-name=blast_MEROPS
#SBATCH --output=logs/MEROPS_split.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=1:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

##this script splits the peptidases from the inhibitors in the MEROPS scan database

#define variables
MEROPS_DB=/hpc/group/vilgalyslab/lal76/databases

#set innput file
input_file="${MEROPS_DB}/merops_scan.lib"

#set output files
peptidase_output="${MEROPS_DB}/merops_scan_peptidase_only.lib"
inhibitor_output="${MEROPS_DB}/merops_scan_inhibitor_only.lib"

#initialize counters
peptidase_count=0
inhibitor_count=0

#set a flag to track if we're in a sequence
in_sequence=false

#loop over the input file
while read -r line; do
    #check that the line starts with '>'
    if [[ "$line" == \>* ]]; then
        #write the sequence to the appropriate output file
        if $in_sequence; then
            if [[ "$header" == *"peptidase unit"* ]]; then
                echo "$header" >> "$peptidase_output"
                echo "$sequence" >> "$peptidase_output"
                ((peptidase_count++))
            elif [[ "$header" == *"inhibitor unit"* ]]; then
                echo "$header" >> "$inhibitor_output"
                echo "$sequence" >> "$inhibitor_output"
                ((inhibitor_count++))
            fi
        fi

        #reset the sequence and set flag
        sequence=""
        in_sequence=true
        header="$line"
    else
        #append the line to the sequence
        sequence+="$line"
    fi
done < "$input_file"

#write the last sequence to the appropriate output file
if $in_sequence; then
    if [[ "$header" == *"peptidase"* ]]; then
        echo "$header" >> "$peptidase_output"
        echo "$sequence" >> "$peptidase_output"
        ((peptidase_count++))
    elif [[ "$header" == *"inhibitor"* ]]; then
        echo "$header" >> "$inhibitor_output"
        echo "$sequence" >> "$inhibitor_output"
        ((inhibitor_count++))
    fi
fi
