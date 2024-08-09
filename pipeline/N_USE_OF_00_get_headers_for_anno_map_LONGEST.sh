#!/bin/bash
#SBATCH --job-name=19_make_mapping_file_for_anno
#SBATCH --output=logs/19_make_mapping_file_for_anno.%A.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --partition=common

#get headers
awk 'sub(/^>/, "")' Orthologs.Longest.fa > longest_headers.txt

#remove everything after the second "|"
cut -f1,2 -d'|' longest_headers.txt > longest_headers2.txt

#split on the "|"
awk -F"|" '$1=$1' OFS=" " longest_headers2.txt > longest_headers3.txt

#in third column, remove everything before the underscore
awk 'gsub(/^[^_]*_/,"",$3)1' longest_headers3.txt > N_USE_longest_headers_clean.txt

rm longest_headers longest_headers2.txt longest_headers3.txt
