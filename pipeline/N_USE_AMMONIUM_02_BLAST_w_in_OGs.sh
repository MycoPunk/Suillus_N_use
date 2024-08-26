#!/bin/bash
#SBATCH --job-name=blast_refrence_ammonium_genes_by_OG_w_missing
#SBATCH --output=logs/blast_ammonium_by_OG_w_missing.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=3:00:00
#SBATCH --mem=28G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

#define variables
QC_CUTOFF=1e-10
N_RESULTS=10
OUTPUT_LOCATION=/hpc/group/vilgalyslab/lal76/Suillus/N_use

#load modules
#module load Python/3.8.1
source /hpc/group/vilgalyslab/lal76/miniconda3/etc/profile.d/conda.sh
module load NCBI-BLAST

#if [ ! -f $INPUT_LOCATION/reference_ammonium_pathway_seqs.fasta.phr ]; then
#  makeblastdb -dbtype prot -in $INPUT_LOCATION/reference_ammonium_pathway_seqs.fasta
#fi

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
    CPUS=1
fi


#re-run BLASTP with missing truncatied seq that didn't end up in OG 0009316
blastp -query $OUTPUT_LOCATION/OG0009316_plus_truncated_seq.fa -db $OUTPUT_LOCATION/ammonium_3_ref.fa -out OG00096316_AMMONIUM.OUT.BLASTP_HITS_w_truncated_seq -num_threads $CPUS -seg yes -soft_masking true -max_target_seqs $N_RESULTS -evalue $QC_CUTOFF -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
