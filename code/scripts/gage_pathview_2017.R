#2/07/17 
#KEGG pathway mapping with pathview
#script for HPC
# R version 3.1.2

library(gage)
#to install gage, had to open R 3.1.2 and install through biocLite

library(pathview)

require(xlsx)
library(dplyr)

#for RNAseq data: http://bioconductor.org/packages/devel/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

setwd("/bio/tgallagh/RNAseq/data/processed_data/GAGE")
#create kegg dataset for sml
kegg_sml <- kegg.gsets("sml")

######### GAGE / PATHVIEW Analysis following : 
#### http://bioconductor.org/packages/devel/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

#Step 1 - Read counts
P30_P36_counts.kegg_final <- read.table(
  "~/GoogleDrive/RNAseq/EdgeR/data/P30_P36_counts.kegg_final.txt", quote="\"", comment.char="")
colnames(P30_P36_counts.kegg_final) <- c("gene", "P30", "P31", "P32", "P33", "P34",
                                         "P35", "P36", "KO")

#table downloaded from KEGG website with K number and brite hierarchy 
sml_k_convert <- read.delim("~/GoogleDrive/RNAseq/KEGG/data/sml_k_convert.txt", header=FALSE)

#add brite #s to table
joined_P30_P36 <- merge(x=P30_P36_counts.kegg_final, y=sml_k_convert, by.x="KO", by.y="V2" )

#get rid of duplicate brite #s
joined_P30_P36 <-joined_P30_P36[!duplicated(joined_P30_P36[,c('V1')]),]
counts_all <- joined_P30_P36[,3:9]


#Step 2 - preprocess counts
#need to normalize # counts by library size
libsizes = colSums(counts_all)

#library sizes is total counts that align to genes
#take log10
log(libsizes)

#average of log libsizes (per million)
mean(log10(libsizes))

size.factor=libsizes/exp(mean(log(libsizes)))

cnts.norm=t(t(counts_all)/size.factor)
cnts.norm=log2(cnts.norm+8)
rownames(cnts.norm) <- joined_P30_P36$V1


cnts_57 <- cnts.norm[,1:5]



#Step 3 = gage function 
ref.idx=3:5
samp.idx=1:2
cnts.kegg.p <- gage(cnts_57, gsets = kegg_sml$kg.sets, ref = ref.idx,
                    samp = samp.idx, compare ="unpaired")
cnts.kegg.p
head(cnts.kegg.p$greater)
heatmap<-sigGeneSet(cnts.kegg.p, heatmap=TRUE)
#Step 4 - pathview 
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
library(pathview)
pv.out.list <- sapply(path.ids2, function(pid) pathview( 
  gene.data = cnts.d, pathway.id = pid,
  species = "sml"))
pv.out <- pathview(gene.data = cnts.d, pathway.id = "sml00270", 
                   species = "sml", gene.idtype="kegg")
################
#repeat for pH 7 versus pH 9 
cnts_79 <-  cnts.norm[,3:7]
#Step 3 = gage function 
ref.idx=1:3
samp.idx=4:5
cnts.kegg.p <- gage(cnts_79, gsets = kegg_sml$kg.sets, ref = ref.idx,
                    samp = samp.idx, compare ="unpaired")
cnts.kegg.p
head(cnts.kegg.p$greater)
heatmap<-sigGeneSet(cnts.kegg.p, heatmap=TRUE)


cnts_59<-  cnts.norm[,-(3:5)]
#Step 3 = gage function 
ref.idx=1:2
samp.idx=3:4
cnts.kegg.p <- gage(cnts_59, gsets = kegg_sml$kg.sets, ref = ref.idx,
                    samp = samp.idx, compare ="unpaired")
cnts.kegg.p
head(cnts.kegg.p$greater)
heatmap<-sigGeneSet(cnts.kegg.p, heatmap=TRUE)

##### RUNS using logFC from EDGE OUTPUT 
### Resulted in no sig upreg pathways...
de5_7_filter.KEGG.final <- read.csv(
  "~/GoogleDrive/RNAseq/EdgeR/data/de5_7_filter.KEGG.final.txt", sep="")
sml_k_convert <- read.delim("~/GoogleDrive/RNAseq/KEGG/data/sml_k_convert.txt", header=FALSE)
joined_5_7 <- merge(x=de5_7_filter.KEGG.final, y=sml_k_convert, by.x="KEGG", by.y="V2" )
joined_5_7<-joined_5_7[!duplicated(joined_5_7[,c('V1')]),]
foldchanges <- joined_5_7$logFC
names(foldchanges) <- joined_5_7$V1

#The gage() function requires a named vector of fold changes, 
#where the names of the values are the gene IDs.
#https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/ 
keggres <- gage(foldchanges, 
                gsets = kegg_sml$kg.sets, same.dir=TRUE)
sigGeneSet(keggres, cutoff =0.9, heatmap=TRUE)
head(keggres$greater)

### 5 vs 9
de5_9_filter.KEGG.final <- read.csv(
  "~/GoogleDrive/RNAseq/EdgeR/data/de5_9_filter.KEGG.final.txt", sep="")
joined_5_9 <- merge(x=de5_9_filter.KEGG.final, y=sml_k_convert, by.x="KEGG", by.y="V2" )
joined_5_9<-joined_5_9[!duplicated(joined_5_9[,c('V1')]),]
foldchanges <- joined_5_9$logFC
names(foldchanges) <- joined_5_9$V1

#The gage() function requires a named vector of fold changes, 
#where the names of the values are the gene IDs.
#https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/ 
keggres <- gage(foldchanges, 
                gsets = kegg_sml$kg.sets, same.dir=TRUE)
sigGeneSet(keggres, cutoff =0.9, heatmap=TRUE)

head(keggres$greater)
