#!/bin/bash
#$ -N bowtie2 
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
ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome

## already built index in command line
#bowtie2-build $ref/Sm_genome.final.scaffolds_CUT.fasta $ref/Sm_index

DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/bowtie2/
bowtie2 --very-fast -I 0 --no-unal -f -x $ref/Sm_index -U $input  -S $DEST$input\.sam
samtools view -bS $DEST$input\.sam | samtools sort > $DEST$input\.sorted.bam
