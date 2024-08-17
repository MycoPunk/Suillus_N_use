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

mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/OG0002448_w_ref.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/nitrate_reductase_alignment_w_ref
mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/OG0002188_w_ref.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/nitrate_transporter_alignment_w_ref
mafft --maxiterate 1000 --localpair --thread 1 /hpc/group/vilgalyslab/lal76/Suillus/N_use/OG0003388_w_ref.fa > /hpc/group/vilgalyslab/lal76/Suillus/N_use/nitrite_reductase_alignment_w_ref
