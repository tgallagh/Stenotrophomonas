#!/bin/bash
#$ -N bowtie2fasta
#$ -q pub64,bio
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-28

#### script for doing bowtie2 alignments of 454 metagenome and metatranscriptome reads against our Sm genome
#load modules

module load bowtie2/2.2.7
module load samtools/1.3

cd  /bio/tgallagh/Stenotrophomonas/data/processed/CFdata

input=$(head -n $SGE_TASK_ID filenames.text | tail -n 1)
ref=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/bowtie2/Smindex/Smindex

## already built index in command line

DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/bowtie2/
bowtie2 --very-fast -I 0 --no-unal -f -x $ref -U $input  -S $DEST$input\.fastafile.sam
samtools view -bS $DEST$input\.fastafile.sam | samtools sort > $DEST$input\.fastafile.sorted.bam
