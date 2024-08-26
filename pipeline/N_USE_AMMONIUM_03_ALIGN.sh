#!/bin/bash
#SBATCH --job-name=align_w_MAFFT
#SBATCH --output=logs/MAFFT.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=1:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --partition=common


#module load MAFFT/v7.543
module unload MAFFT
module load MAFFT/7.475-rhel8

#mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/OG0009316_w_ref_QC.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/ammonium_OG0009316_alignment_w_ref.aln
#mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/OG0006295_w_ref_QC.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/ammonium_OG0006295_alignment_w_ref.aln
#mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/OG0006824_w_ref_QC.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/ammonium_OG0006824_alignment_w_ref.aln

mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/ammonium_all_3_ref_seqs.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/ammonium_all_3_ref_alignment.aln
