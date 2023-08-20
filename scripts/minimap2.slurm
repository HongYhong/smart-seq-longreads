#! /bin/bash

#SBATCH -J minimap2
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 80

for i in fastq_hj/*.fastq.gz
do
	minimap2  -t 80 -a -x splice --splice-flank=yes genome/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz  $i |samtools sort -@ 80 > alignments_hj_mm/`basename $i .fastq.gz`.sorted.bam
done
