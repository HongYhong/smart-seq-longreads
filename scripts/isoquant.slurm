#! /bin/bash

#SBATCH -J ISOQUANT
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=100G
#SBATCH --ntasks 50

~/IsoQuant-master/isoquant.py --low_memory -t 50 --output isoquant
--matching_strategy loose \
--sqanti_output --model_construction_strategy sensitive_ont -d nanopore \
--bam_list MFA1912RKSv2_bamlist.txt
--read_group file_name --complete_genedb \
--genedb annotation/GCF_012559485.2_MFA1912RKSv2_genomic.corrected.gtf \
--reference genome/GCF_012559485.2_MFA1912RKSv2_genomic.fna \
--clean_start --force
