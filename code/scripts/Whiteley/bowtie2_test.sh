#!/bin/bash
#$ -N bowtie2htseq
#$ -q pub64,bio
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1

## Script for doing bowtie2 alignments of quality filtered Whiteley reads (using trimmmomatic)
## alignments onto the entire genome 

#load modules
module load bowtie2/2.2.7
module load samtools/1.3

cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/

input=derepderep.Rn_Meta_Invivo_HumanSputum_2016_JLD_HumanDanishSputum_E.trimmed.fastq.fastq
ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome
DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie2/derep/

bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -q -x$ref/Sm_index \
-U  $input \
-S  $DEST$input\.derep.sam

samtools view -bS $DEST$input\.derep.sam | samtools sort > $DEST$input\.derep.sorted.bam

HTSEQ=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/HTSeq/
cd /bio/tgallagh/programs/HTSeq-0.6.1/scripts
samtools view  $DEST$input\.derep.sorted.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ$input\.xls

