#!/bin/sh

# Tara Script to summarize pre-processing of  Sm RNAseq reads
# 5/03/17


#output file
touch /bio/tgallagh/Stenotrophomonas/output/tables/Steno_summary.txt
OUTPUT="/bio/tgallagh/Stenotrophomonas/output/tables/Steno_summary.txt"


echo "#RAW READS SUMMARY" > $OUTPUT
#Determine the number of reads per file
cd /bio/tgallagh/Stenotrophomonas/data/raw/reads/Sm

for i in $(seq 30 36); do 

echo  Raw_Fwd_Read_P$i\ $'\t' "$(grep ^"@" <P$i\_READ1.fastq \
      | wc -l)" >> $OUTPUT;
echo  Raw_Rev_Read_P$i\ $'\t' "$(grep ^"@" <P$i\_READ2.fastq \
      | wc -l)" >> $OUTPUT
done

#Determine number of reads that pass trimmomatic QF
#change directory
cd /bio/tgallagh/Stenotrophomonas/data/processed/trimmed

echo "#TRIMMOMATIC OUTPUT SUMMARY" >> $OUTPUT

echo "#READS WHERE BOTH MATES PASSES QF" >> $OUTPUT
for i in $(seq 30 36); do
echo  Trimm_Fwd_PairedRead_P$i\ $'\t' "$(grep ^"@" <P$i\_READ1_paired2.fastq \
      | wc -l)" >> $OUTPUT;
echo  Trimm_Rev_PairedRead_P$i\ $'\t' "$(grep ^"@" <P$i\_READ2_paired2.fastq \
      | wc -l)" >> $OUTPUT
done

echo "#READS WHERE 1 MATE PASSES QF" >> $OUTPUT

for i in $(seq 30 36); do
echo  Trimm_Fwd_Read_P$i\ $'\t' "$(grep ^"@" <P$i\_READ1_unpaired2.fastq \
      | wc -l)" >> $OUTPUT;
echo  Trimm_Rev_Read_P$i\ $'\t' "$(grep ^"@" <P$i\_READ2_unpaired2.fastq \
      | wc -l)" >> $OUTPUT
done 

#Determine number of reads that are overlapping
echo "#OVERLAPPING PE READS" >> $OUTPUT

for i in $(seq 30 36); do
echo  Overlapping_reads_P$i\ $'\t' "$(grep ^"@" <P$i\_READ1_READ2_merged.fastq.assembled.fastq\
| wc -l)" >> $OUTPUT;
done

#Determine number of reads with rRNA contaminants
echo "#FINAL CLEANED READS AFTER RRNA REMOVAL" >> $OUTPUT

for i in $(seq 30 36); do
echo SE_Final_P$i\ $'\t' "$(grep ^">" <P$i\_READ1_READ2_cat_unpaired.deconseq.clean.fasta \
| wc -l)" >> $OUTPUT
echo  Forward_Final_P$i\ $'\t' "$(grep ^">" <P$i\_READ1_READ2_merged.fastq.unassembled.forward.deconseq.clean.paired.fasta \
| wc -l)" >> $OUTPUT
echo Reverse_Final_P$i\ $'\t' "$(grep ^">" <P$i\_READ1_READ2_merged.fastq.unassembled.reverse.deconseq.clean.paired.fasta \ 
| wc -l)" >> $OUTPUT
done


### Determine number of bowtie2 alignments
module load samtools
 
cd /bio/tgallagh/Stenotrophomonas/data/processed/bowtie2/raw/
echo "PERCENT RAW READS THAT ALIGNED" >> $OUTPUT 
for i in $(seq 30 36); do
samtools flagstat P$i\_raw.sorted.bam >  P$i\.raw.stats.temp
echo P$i\_raw_alignment $(cat P$i\.raw.stats.temp | head -5 | tail -1) >> $OUTPUT
done

cd /bio/tgallagh/Stenotrophomonas/data/processed/bowtie2/
echo "PERCENT PROCESSED PAIRED READS THAT ALIGNED" >> $OUTPUT

for i in $(seq 30 36); do
samtools flagstat P$i\_both.sorted.bam > P$i\_paired.stats.temp
echo Final_Paired_Align_P$i\ $'\t' $(cat P$i\_paired.stats.temp | head -5 | tail -1) >> $OUTPUT
done

echo "PERCENT PROCESSED OVERLAPPED READS (SE) THAT ALIGNED" >> $OUTPUT
for i in $(seq 30 36); do
samtools flagstat P$i\_singles.sorted.bam > P$i\_singles.stats.temp
echo Final_Overlapped_Align_P$i\ $'\t' $(cat P$i\_singles.stats.temp | head -5 | tail -1) >> $OUTPUT
done

rm -f *stats.temp




