#!/bin/bash
#SBATCH --job-name=blast_Lbicolor_ammonium_genes
#SBATCH --output=logs/blast_ammonium.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=3:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

#define variables
OUTEXT=blasttab
EXT=aa.fasta
INFILE=Orthologs.Longest.fa
QC_CUTOFF=1e-10
N_RESULTS=10
INPUT_LOCATION=/hpc/group/vilgalyslab/lal76/Suillus/N_use


#load modules
#module load Python/3.8.1
source /hpc/group/vilgalyslab/lal76/miniconda3/etc/profile.d/conda.sh
module load NCBI-BLAST

if [ ! -f $INPUT_LOCATION/ammonium_3_ref.fa.phr ]; then
  makeblastdb -dbtype prot -in $INPUT_LOCATION/ammonium_3_ref.fa
fi

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
    CPUS=1
fi


#w/query length to calculate % coverage
blastp -query $INFILE -db $INPUT_LOCATION/ammonium_3_ref.fa -out AMMONIUM.3_ref_OUT.BLASTP_HITS.aa.fasta -num_threads $CPUS -seg yes -soft_masking true -max_target_seqs $N_RESULTS -evalue $QC_CUTOFF -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
