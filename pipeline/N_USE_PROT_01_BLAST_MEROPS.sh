#!/bin/bash
#SBATCH --job-name=blast_MEROPS
#SBATCH --output=logs/blast_MEROPS.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=3:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

#this script uses BLASTP to search MEROPS proteases against the longest representative of each OG

#define variables
OUTEXT=blasttab
EXT=aa.fasta
INFILE=Orthologs.Longest.fa
QC_CUTOFF=1e-10
N_RESULTS=10
MEROPS_DB=/hpc/group/vilgalyslab/lal76/databases


#load modules
#module load Python/3.8.1
source /hpc/group/vilgalyslab/lal76/miniconda3/etc/profile.d/conda.sh
module load NCBI-BLAST

if [ ! -f $MEROPS_DB/merops_scan_peptidase_only.lib.phr ]; then
  makeblastdb -dbtype prot -in $MEROPS_DB/merops_scan_peptidase_only.lib
fi

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
    CPUS=1
fi

#run blastp
blastp -query $INFILE -db $MEROPS_DB/merops_scan_peptidase_only.lib -out MEROPS.OUT.PEPTIDASE.aa.fasta -num_threads $CPUS -seg yes -soft_masking true -max_target_seqs $N_RESULTS -evalue $QC_CUTOFF -use_sw_tback -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
