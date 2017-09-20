#!/bin/bash
#$ -N bowtie2_Whiteley
#$ -q pub64,bio,free64
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7

#### script for doing bowtie2 alignments of Whiteley metatranscriptome reads against our Sm genome and against Sm DGE

#load modules
module load bowtie2/2.2.7
module load samtools/1.3

cd  /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/fastq

input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome
ref2=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/454data/Smindex

## already built indexes for both genome and fasta of DGE in command line

DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie2/
bowtie2 --very-fast -I 0 --no-unal -q -x $ref/Sm_index -U $input  -S $DEST$input\.sam
samtools view -bS $DEST$input\.sam | samtools sort > $DEST$input\.sorted.bam

DEST2=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie2/
bowtie2 --very-fast -I 0 --no-unal -q -x $ref2/Smindex -U $input  -S $DEST2$input\.fastafile.sam
samtools view -bS $DEST2$input\.fastafile.sam | samtools sort > $DEST2$input\.fastafile.sorted.bam

