#!/bin/bash
#$ -N Sm_RNAseq_pipeline
#$ -q pub64,bio
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-7

#### Pipeline for quality filtering and aligning Sm RNAseq data

### 

#load modules
#have to purge before using QIIME
#make sure no other python loaded

module purge
module load enthought_python/7.3.2
module load PEAR/20150409
module load bowtie2/2.2.7
module load samtools/1.3
module load trimmomatic/0.35
module load perl
module load fastx_toolkit/0.0.14
module load jje/jjeutils/0.1a
module load bamtools/2.3.0
module load bedtools/2.25.0

### Name task IDs by prefix of file names 
cd  /bio/tgallagh/Stenotrophomonas/data/raw/reads/Sm


#Prepare the prefix text file in interactive mode:
#for x in *READ1.fastq; do echo $x >> prefix.txt; done
#sed -i "s/"READ1.fastq"//g"  prefix.txt

input=$(head -n $SGE_TASK_ID prefix.txt | tail -n 1)

ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome


##################################
##### Look at Raw Reads alignments
##################################
### Sloppy bowtie2 alignments of raw reads
#bowtie2-build $ref/FLR01_A5.final.scaffolds.cut.fasta $ref/P1index

RAWDEST=/bio/tgallagh/Stenotrophomonas/data/processed/bowtie2/raw/

## build index in command line
#bowtie2-build $ref/Sm_genome.final.scaffolds_CUT.fasta $ref/Sm_index

bowtie2 --very-fast -x $ref/Sm_index -1 $input\READ1.fastq -2 $input\READ2.fastq  -S $RAWDEST$input\raw.sam
samtools view -bS $RAWDEST$input\raw.sam | samtools sort > $RAWDEST$input\raw.sorted.bam

#############################
##### Clean up the raw reads
#############################

#### A lot of quality filtering! 

### High adaptor content from fastqc
# Use Trimmomatic to remove adaptors and quality filter

ADAPTORS=/bio/tgallagh/Stenotrophomonas/data/raw/adaptors/

#mkdir /bio/tgallagh/Stenotrophomonas/data/processed/trimmed
TRIMMED=/bio/tgallagh/Stenotrophomonas/data/processed/trimmed/

java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE -phred33 $input\READ1.fastq $input\READ2.fastq $TRIMMED$input\READ1_paired2.fastq $TRIMMED$input\READ1_unpaired2.fastq $TRIMMED$input\READ2_paired2.fastq $TRIMMED$input\READ2_unpaired2.fastq ILLUMINACLIP:$ADAPTORS$input\CustomPrimers.sh:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:4:20 MINLEN:50

### 250 PE reads so probably a lot of overlap
# use PEAR to merge overlaps from trimmomatic output
pear -f $TRIMMED$input\READ1_paired2.fastq -r $TRIMMED$input\READ2_paired2.fastq -o $TRIMMED$input\READ1_READ2_merged.fastq

#remove .discarded files
rm -f $TRIMMED/$inputREAD1_READ2_merged.fastq.discarded.fastq

### prepare reads for next set of steps
# convert PEAR output fastq to fasta
fastq_to_fasta -n -i $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.fastq -o $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.fasta
fastq_to_fasta -n -i $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.fastq -o $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.fasta
#cat single mates from the TRIMMOMATIC STEP
cat $TRIMMED$input\READ1_READ2_merged.fastq.assembled.fastq $TRIMMED$input\READ1_unpaired2.fastq  $TRIMMED$input\READ2_unpaired2.fastq  > $TRIMMED$input\READ1_READ2_cat_unpaired2.fastq
fastq_to_fasta -n -i $TRIMMED$input\READ1_READ2_cat_unpaired2.fastq -o $TRIMMED$input\READ1_READ2_cat_unpaired2.fasta

### Shouldn't be a lot of rRNA since we removed with kit
# just incase, use deconseq to remove rRNA contaminants
# already made database consisting of P1 rRNA annotated by RAST
cd /bio/tgallagh/programs/deconseq
#have to change directory to get deconseq to run
perl /bio/tgallagh/programs/deconseq/deconseq.pl -dbs SmrRNA  -f $TRIMMED$input\READ1_READ2_cat_unpaired2.fasta -i 95 -c 90 -o  $TRIMMED$input\READ1_READ2_cat_unpaired2.fasta.deconseq
perl /bio/tgallagh/programs/deconseq/deconseq.pl -dbs SmrRNA -f $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.fasta  -i 95 -c 90 -o $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.fasta.deconseq
perl /bio/tgallagh/programs/deconseq/deconseq.pl -dbs SmrRNA -f $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.fasta -i 95 -c 90 -o  $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.fasta.deconseq

### change file names
cd $TRIMMED$input\READ1_READ2_cat_unpaired2.fasta.deconseq
cp *clean.fa $TRIMMED$input\READ1_READ2_cat_unpaired.deconseq.clean.fasta
cp *cont.fa $TRIMMED$input\READ1_READ2_cat_unpaired.deconseq.contam.fasta

cd $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.fasta.deconseq
cp *clean.fa $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.clean.fasta
cp *cont.fa $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.contam.fasta

cd $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.fasta.deconseq
cp *clean.fa $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.clean.fasta
cp *cont.fa $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.contam.fasta

### Since we filtered the forward and reverse reads seperately, need to resync these files
bioawk -c fastx '{print $name } ' \
< $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.contam.fasta \
> $TRIMMED$input\reverse_contam_list
bioawk -c fastx '{print $name } ' \
< $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.contam.fasta \
> $TRIMMED$input\forward_contam_list
#make list of rRNA contam reads from both forward and reverse reads files
cat $TRIMMED$input\reverse_contam_list $TRIMMED$input\forward_contam_list > $TRIMMED$input\both_contam_list
#unique command to make sure no double contaminants
cat $TRIMMED$input\both_contam_list | sort | uniq > $TRIMMED$input\both_contam_list_uniq
#use QIIME python script to filter out unwanted sequences
#should have equal # of reads and (same read identifier) in forward and reverse reads files now
filter_fasta.py -f $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.clean.fasta \
-o $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.clean.paired.fasta -s $TRIMMED$input\both_contam_list_uniq -n 
filter_fasta.py -f $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.clean.fasta \
-o $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.clean.paired.fasta -s $TRIMMED$input\both_contam_list_uniq -n

#############################
#### Align reads with bowtie2
#############################
# the overlapping reads are input as SE reads
ALIGN=/bio/tgallagh/Stenotrophomonas/data/processed/bowtie2/
bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -f -x$ref/Sm_index \
-U  $TRIMMED$input\READ1_READ2_cat_unpaired.deconseq.clean.fasta \
-S  $ALIGN$input\singles.sam
# the filtered, resync PE reads are input as PE reads
bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -f -x$ref/Sm_index \
-1  $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.forward.deconseq.clean.paired.fasta \
-2  $TRIMMED$input\READ1_READ2_merged.fastq.unassembled.reverse.deconseq.clean.paired.fasta \
-S  $ALIGN$input\both.sam

#use Samtools to convert to bam and sort the alignment files
samtools view -bS $ALIGN$input\singles.sam | samtools sort > $ALIGN$input\singles.sorted.bam
samtools view -bS $ALIGN$input\both.sam | samtools sort > $ALIGN$input\both.sorted.bam

##########################################
# Option 1 - add read groups for use later
##########################################
#java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=$ALIGN$input\singles.sorted.bam O=$ALIGN$input\singles.sorted.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID={$input} RGSM={$input} VALIDATION_STRINGENCY=LENIENT
#java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=$ALIGN$input\both.sorted.bam O=$ALIGN$input\both.sorted.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID={$input} RGSM={$input} VALIDATION_STRINGENCY=LENIENT
############################################################
# Option 2 - use htseq to get read counts to annotated genes
############################################################
#mkdir/bio/tgallagh/EE283/EE283_final/data/processed/htseq
HTSEQ=/bio/tgallagh/Stenotrophomonas/data/processed/HTSeq/
cd /bio/tgallagh/programs/HTSeq-0.6.1/scripts
samtools view $ALIGN$input\singles.sorted.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ$input\single.xls

samtools view $ALIGN$input\both.sorted.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ$input\both.xls

#take sum of PE and SE files 
paste $HTSEQ$input\both.xls $HTSEQ$input\singles.xls | cut -f-1,2,4 | awk '{print $2 + $3}'>  $HTSEQ$input\.temp
paste $HTSEQ$input\both.xls $HTSEQ$input\.temp | cut -f-1,3  > $HTSEQ$input\.2.temp

##########################################################
# Option 3 prepare coverage files to look at coverage in R 
#########################################################
bedtools genomecov -ibam $ALIGN$input\singles.sorted.bam -g /bio/tgallagh/Stenotrophomonas/data/processed/genome/Sm_cut.genome -d > $ALIGN$input\singles.sorted.coverage.txt
bedtools genomecov -ibam $ALIGN$input\both.sorted.bam -g /bio/tgallagh/Stenotrophomonas/data/processed/genome/Sm_cut.genome -d > $ALIGN$input\both.sorted.coverage.txt
#get sum of two files
paste $ALIGN$input\singles.sorted.coverage.txt $ALIGN$input\both.sorted.coverage.txt | cut -f-1,2,3,6 | awk '{print $3+$4}' >  $ALIGN$input\temp.sum
mkdir /bio/tgallagh/Stenotrophomonas/data/processed/bowtie2/coverage
COV=/bio/tgallagh/Stenotrophomonas/data/processed/bowtie2/coverage/
paste $ALIGN$input\singles.sorted.coverage.txt $ALIGN$input\temp.sum | cut -f-1,2,4 > $COV$input\final.coverage.txt





