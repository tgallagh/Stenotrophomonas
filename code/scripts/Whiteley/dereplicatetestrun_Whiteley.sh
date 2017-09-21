#!/bin/bash
#$ -N WhiteleyDerepTestRun
#$ -q free64,pub64,bio
#$ -m e

#### Dereplicate with Whiteley Sputum E sample (quality filtered)
#### 09/21/2017
#### Look at effect of the number of duplicates on downstream DGE

# path to quality-filtered, trimmed SE reads for sample E
INPUT=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/trimmed/Rn_Meta_Invivo_HumanSputum_2016_JLD_HumanDanishSputum_E.trimmed.fastq
#get rid of ALL exact duplicates
module load prinseq-lite/0.20.4 
DEREPTEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/derep/testrun/
prinseq-lite.pl -fastq $INPUT -derep 1 -out_good $DEREPTEST\derep1.$INPUT -out_bad null -out_format 3
#prinseq report
REPORT=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/derep/testrun/prinseqreports/
prinseq-lite.pl -fastq $DEREPTEST\derep1.$INPUT -graph_data $REPORT\derep1.$INPUT\.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn
prinseq-graphs.pl -i $REPORT\derep1.$INPUT\.gd -html_all -o $REPORT
#bowtie2 alignment
module load bowtie2/2.2.7
module load samtools/1.3
ref=/bio/tgallagh/Stenotrophomonas/data/processed/genome
DEST=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/derep/testrun/alignments/

bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -q -x$ref/Sm_index \
-U  $DEREPTEST\derep1.$INPUT \
-S  $DEST$input\derep1.$INPUT\.sam
samtools view -bS $DEST$input\derep1.$INPUT\.sam | samtools sort >$DEST$input\derep1.$INPUT\.bam

# htseq count
HTSEQ=/bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/derep/testrun/alignments/
cd /bio/tgallagh/programs/HTSeq-0.6.1/scripts
samtools view $DEST$input\derep1.$INPUT\.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ\.derep1.xls

### repeat for 3
prinseq-lite.pl -fastq $INPUT -derep 3 -out_good $DEREPTEST\derep3.$input -out_bad null -out_format 3
prinseq-lite.pl -fastq $DEREPTEST\derep3.$INPUT -graph_data $REPORT\derep3.$INPUT\.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn
prinseq-graphs.pl -i $REPORT\derep3.$INPUT\.gd -html_all -o $REPORT

bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -q -x$ref/Sm_index \
-U  $DEREPTEST\derep3.$INPUT \
-S  $DEST$input\derep3.$INPUT\.sam
samtools view -bS $DEST$input\derep3.$INPUT\.sam | samtools sort >$DEST$input\derep3.$INPUT\.bam

cd /bio/tgallagh/programs/HTSeq-0.6.1/scripts
samtools view $DEST$input\derep3.$INPUT\.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ\.derep3.xls

## repeat for 10
prinseq-lite.pl -fastq $INPUT -derep 10 -out_good $DEREPTEST\derep10.$input -out_bad null -out_format 3
prinseq-lite.pl -fastq $DEREPTEST\derep10.$INPUT -graph_data $REPORT\derep10.$INPUT\.gd -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn
prinseq-graphs.pl -i $REPORT\derep10.$INPUT\.gd -html_all -o $REPORT

bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I 0 --no-unal -q -x$ref/Sm_index \
-U  $DEREPTEST\derep10.$INPUT \
-S  $DEST$input\derep10.$INPUT\.sam
samtools view -bS $DEST$input\derep10.$INPUT\.sam | samtools sort >$DEST$input\derep10.$INPUT\.bam

cd /bio/tgallagh/programs/HTSeq-0.6.1/scripts
samtools view $DEST$input\derep10.$INPUT\.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ\.derep10.xls

## without any duplicate removal
# already generated alignment file for trimmomatic output
cd /bio/tgallagh/programs/HTSeq-0.6.1/scripts
samtools view /bio/tgallagh/Stenotrophomonas/data/processed/CFdata/Whiteley/bowtie2/qualityfiltered/Rn_Meta_Invivo_HumanSputum_2016_JLD_HumanDanishSputum_E.fastq.trimmed.trimmed.sorted.bam | ./htseq-count -s no -t CDS -i ID - /bio/tgallagh/Stenotrophomonas/data/processed/genome/RAST/6666666.230262.gff > $HTSEQ\.noderep.xls


