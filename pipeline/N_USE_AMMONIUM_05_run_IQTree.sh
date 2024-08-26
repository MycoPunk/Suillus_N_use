#!/usr/bin/bash
#SBATCH -p common -n 4 -N 1 --mem 16gb --out iqtree_gene.%A.log

module load IQ-TREE/1.6.12

iqtree -s ammonium_all_3_ref_alignment.aln.clipkit -alrt 1000 -bb 1000 -m MFP
