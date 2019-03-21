#!/bin/bash
#$ -N bowtie2Database
#$ -q pub64,bio,free64
#$ -pe make 8 
#$ -R y
#$ -t 1-13

#script for doing bowtie2 alignments on CF genome database, check sensitivity of bowtie2 parameters for Steno metatranscriptomes

#load modules
module load bowtie/1.0.0
module load samtools/1.3

# cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed

#cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdatabase/genomes

#input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
#ref=/bio/tgallagh/Stenotrophomonas/data/processed/CFdatabase/genomes/
#DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/database/
#bowtie2-build $ref$input $input
#bowtie-build $ref$input $input

#bowtie -v 0 -l 75 -y --best -q $ref$input /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie1/sputum.mapped.pooled.fastq \
-S $DEST$input\.sam

### HTSEQ COUNT
module load samtools/1.3

cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/database
HTSEQ=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/database/

input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
samtools view $input | /bio/tgallagh/programs/HTSeq-0.6.1/scripts/htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ$input\.xls


