#! /bin/bash

#SBATCH -J salmon_index
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks 10

salmon index -t genome/Homo_sapiens.GRCh38.cdna.all.fa.gz -i Homo_GRCh38_salmon_index/Homo_GRCh38_salmon_index







