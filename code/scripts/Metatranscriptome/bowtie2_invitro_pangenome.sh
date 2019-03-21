#!/bin/bash
#$ -N P33
#$ -q pub64,bio,free64
#$ -pe make 8 
#$ -R y
#$ -t 1-7

#load modules
module load bowtie2/2.2.3 
module load samtools/1.3


## bowtie2-build POOLED.CF.NO.STENO.fna POOLED.CF
### align the in vitro reads against our CF database
#input=$(head -n $SGE_TASK_ID /bio/tgallagh/Stenotrophomonas/data/raw/reads/Sm/prefix.txt | tail -n 1)
#ref=/bio/tgallagh/databases/CFdatabase_noSteno/
#DEST=/bio/tgallagh/Stenotrophomonas/data/processed/invitro_database_align/
#TRIMMED=/bio/tgallagh/Stenotrophomonas/data/processed/trimmed/
overlapping reads from quality filtering are input as SE reads


#bowtie2 --very-sensitive --un /bio/tgallagh/Stenotrophomonas/data/processed/invitro_database_unaligned/$input\singles.fasta -x $ref\POOLED.CF -f -U $TRIMMED$input\READ1_READ2_cat_unpaired.deconseq.clean.fasta  -S $DEST$input\singles.sam
## both mates there == > Paired end alignement in bowtie2 

#bowtie2 --maxins 2000 --very-sensitive --un-conc /bio/tgallagh/Stenotrophomonas/data/processed/invitro_database_unaligned/$input\paired.fasta -x $ref\POOLED.CF -f -1  $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.clean.paired.fasta \
#-2  $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.clean.paired.fasta  \
#-S $DEST$input\paired.sam



### align the reads that didn't align to database to the pan genome
cd  /bio/tgallagh/Stenotrophomonas/data/processed/invitro_database_unaligned/
ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome/phylogenetic/tara_fasta/fixed_fasta/POOLED/POOLED.STENO.CDS
DEST=/bio/tgallagh/Stenotrophomonas/data/processed/invitro_align_final/

bowtie2 -a --score-min C,0,0 -x $ref -f P33_singles.fasta -S $DEST$\P33_singles.sam

bowtie2 -a --maxins 1000 --score-min C,0,0 -x $ref -f -1 P33_paired.1.fasta -2 P33_paired.2.fasta -S $DEST\P33_paired.sam

