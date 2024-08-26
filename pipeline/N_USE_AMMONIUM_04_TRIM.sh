#!/bin/bash
#SBATCH --job-name=trim_alignment
#SBATCH --output=logs/trim_alignment.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=3:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

module load Python

source .venv/bin/activate

clipkit ammonium_all_3_ref_alignment.aln
