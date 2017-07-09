
#This folder contains HTseq output aligning to RAST annotation

# for each sample, need to take the sum of "single" and "both"
	# single = SE reads
	# both = both PEs used in alignment
for i in $(seq 30 36);
do 
paste P$i\_single.xls P$i\_both.xls | awk '{print $1, ($2 + $4) }' > P$i\_sum.xls
done

#combine all sample HTSEQ input into one table
for i in $(seq 30 36);
do
cat P$i\_sum.xls | cut -d ' ' -f 2  > P$i\_sum.temp
done

paste P30_sum.xls P31_sum.temp P32_sum.temp P33_sum.temp P34_sum.temp P35_sum.temp P36_sum.temp > steno_htseq_count.temp

sed -i "1igene P30 P31 P32 P33 P34 P35 P36" steno_htseq_count.temp | column -t > stenohtseq_final.xls 
column -t steno_htseq_count.temp > steno_htseq.xls



#annotate counts with KEGG number
join -1 1 -2 1 -o 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 1.2 <(sort -k 1b,1 /bio/tgallagh/RNAseq/data/annotations/KEGG_Sm_CUT_sml.ko) <(sort -k 1b,1 P30_P36_counts_final.txt) | column -t > P30_P36_counts.kegg.txt 
grep K P30_P36_counts.kegg.txt > P30_P36_counts.kegg_final.txt

`
:
