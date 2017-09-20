#!/bin/bash
#$ -N bowtie2WhiteleyQF
#$ -q pub64,bio,free64
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7

## Script for doing bowtie2 alignments of quality filtered Whiteley reads (using trimmmomatic)
## alignments onto the entire genome 

#load modules
module load bowtie2/2.2.7
module load samtools/1.3

cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed

input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome
DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie2/qualityfiltered/

bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -q -x$ref/Sm_index \
-U  $input \
-S  $DEST$input\.trimmed.sam

samtools view -bS $DEST$input\.trimmed.sam | samtools sort > $DEST$input\.trimmed.sorted.bam

