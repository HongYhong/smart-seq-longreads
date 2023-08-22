#! /bin/bash

#SBATCH -J smartseq-salmon
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks 50
#SBATCH --array=0-12

readarray -t fqs1 < illumina/first_fq.txt
readarray -t fqs2 < illumina/second_fq.txt

index_file="Homo_GRCh38_salmon_index/Homo_GRCh38_salmon_index"

fq1=${fqs1[$SLURM_ARRAY_TASK_ID]}
fq2=${fqs2[$SLURM_ARRAY_TASK_ID]}

sample=`basename $fq1 .fq.gz`

salmon_cmd="salmon quant -i $index_file -l A \
         -1 $fq1 \
         -2 $fq2 \
         -p 50 -o Homo_GRCh38_salmon_quant/${sample}_quant"

eval ${salmon_cmd}

