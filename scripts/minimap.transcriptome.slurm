#! /bin/bash

#SBATCH -J minimap2_transcirptome
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 80

for i in fastq_hj/*.fastq.gz
do
	minimap2  -t 80 -ax map -N 100 -p 0.99  genome/Homo_sapiens.GRCh38.cdna.all.fa.gz $i |samtools sort -@ 80 > alignment_hj_transcriptome/`basename $i .fastq.gz`.sorted.bam
done
