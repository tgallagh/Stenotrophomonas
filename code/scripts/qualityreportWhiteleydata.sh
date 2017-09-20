#!/bin/bash
#$ -N qualitytrimmedset1Whitley
#$ -q pub64,bio,free64
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7



# script for prinseq-lite quality report of raw Whiteley Mtr reads
module load prinseq-lite/0.20.4 

cd  /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/fastq

# For raw reads
#input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
#DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/454data/prinseqreports/rawreads/
#prinseq-lite.pl -fastq $input -graph_data $DEST$input\.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn

# preinseq report of quality filtered reads (set 1) 
FILENAMES=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed/filenames.txt
input=$(head -n $SGE_TASK_ID $FILENAMES | tail -n 1)
cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed/

DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/454data/prinseqreports/qualityfilteredset1
prinseq-lite.pl -fastq $input -graph_data $DEST$input\.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn


