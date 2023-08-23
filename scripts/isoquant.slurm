#! /bin/bash

#SBATCH -J ISOQUANT
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 50

/cluster/home/hongyanhong/IsoQuant-master/isoquant.py --low_memory -t 50 --output control_data/isoquant \
--matching_strategy loose \
--sqanti_output --model_construction_strategy sensitive_ont -d nanopore \
--bam_list control_data/isoquant/bamlist.txt \
--read_group file_name --complete_genedb \
--genedb genome/Homo_sapiens.GRCh38.109.chr.gtf \
--reference genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--clean_start --force
