#!/bin/bash
#$ -N WhiteleyDerep
#$ -q free64,pub64,bio
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7

#### quality filter Whiteley reads with trimmomatic

module load prinseq-lite/0.20.4 
cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed/
input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)

TRIMMED=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/derep/
prinseq-lite.pl -fastq $input -derep 1 -derep_min 2  -out_good $TRIMMED\derep.$input -out_bad null -out_format 3 


