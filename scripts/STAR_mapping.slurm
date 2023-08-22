#! /bin/bash

#SBATCH -J STAR_mapping
#SBATCH -p q1,q2,q3,fat01,gpu1
#SBATCH -N 1
#SBATCH --output=log.%A.out
#SBATCH --error=log.%A.err
#SBATCH --mem=100G
#SBATCH --ntasks 50
#SBATCH --array=0-12

readarray -t fqs1 < illumina/first_fq.txt
readarray -t fqs2 < illumina/second_fq.txt

fq1=${fqs1[$SLURM_ARRAY_TASK_ID]}
fq2=${fqs2[$SLURM_ARRAY_TASK_ID]}

STAR --runThreadN 50  --genomeDir Homo_GRCh38_STAR --readFilesIn ${fq1} ${fq2} \
 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --readFilesCommand zcat --outFileNamePrefix \
  illumina_aln/`basename $fq1 .fq.gz`.sorted.bam
