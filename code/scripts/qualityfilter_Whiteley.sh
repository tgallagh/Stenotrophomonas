#!/bin/bash
#$ -N WhiteleyQF
#$ -q free64,pub64,bio
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7

#### quality filter Whiteley reads with trimmomatic

module load trimmomatic/0.35

cd  /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/fastq
input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)

TRIMMED=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed/

java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 $input $TRIMMED$input\.trimmed ILLUMINACLIP:/bio/tgallagh/alladaptors.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:4:20 MINLEN:35

