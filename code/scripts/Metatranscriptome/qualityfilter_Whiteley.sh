#!/bin/bash
#$ -N WhitelyFilter
#$ -q free64,pub64,
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7

#### quality filter Whiteley reads with trimmomatic
# module load trimmomatic/0.35

#cd  /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/fastq
#input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
#TRIMMED=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmedfinal/

#java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 $input $TRIMMED$input\.trimmed ILLUMINACLIP:/bio/tgallagh/alladaptors.fa:2:30:10 LEADING:30 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:45

#### Dereplicate duplicates after quality filtering
#cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmedfinal/
#input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
## get rid of ALL exact duplicates
#module load prinseq-lite/0.20.4
#DEREP=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/derep/

#prinseq-lite.pl -fastq $input -derep 1 -out_good $DEREP$input -out_bad null -out_format 3
#prinseq report
#REPORT=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/prinseqreports/derep/
#prinseq-lite.pl -fastq $DEREP$input\.fastq -graph_data $REPORT$input\.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn

## align with bowtie1 
#module load bowtie/1.0.0

#DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie1/

#bowtie -v 0 -l 75 -y --best -q /bio/tgallagh/Stenotrophomonas/data/processed/genome/Sm_index  $DEREP$input\.fastq  \
#-S $DEST$input\.sam

### HTSEQ COUNT
module load samtools/1.3

cd /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie1
HTSEQ=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/htseq/
input=$(head -n $SGE_TASK_ID filenames.txt | tail -n 1)
samtools view $input | /bio/tgallagh/programs/HTSeq-0.6.1/scripts/htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ$input\.xls


