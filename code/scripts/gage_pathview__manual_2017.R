### 8/20/17
#KEGG pathway mapping with pathview
#script for HPC
# R version 3.3.1

library(gage)
#to install gage, had to open R 3.1.2 and install through biocLite

library(pathview)

require(xlsx)
library(dplyr)
#for RNAseq data: http://bioconductor.org/packages/devel/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
setwd("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/")

#create kegg dataset for sml
kegg_sml <- kegg.gsets("sml")




##########################################################
### FISHER EXACT TEST TO LOOK FOR OVER REPRESENTED GROUPS
##########################################################
#helpful slides = http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
# get number of genes from the seed annotation file
SEED <- read.delim("~/GoogleDrive/Stenotrophomonas/output/tables/reformmatedseedannotations.txt")
#this table has duplicated features, so need to get # unique genes 
SEED_nodup <- SEED[!duplicated(SEED$Features),]
number_genes <- nrow(SEED_nodup)

# data frame with all edgeR upreg and downreg genes (pH 5 and pH 9 conditions)
# has RAST SEED categories 
data <- read.xlsx("~/GoogleDrive/Stenotrophomonas/output/tables/edgeRallRAST.xlsx", 1)

#just upreg in pH 5 without uncharacterized genes 
up5 <- subset(data,which=="5up")
up5 <- subset(up5, !Category == "Uncharacterized")

down5 <- subset(data,which=="5down")
down5 <- subset(down5, !Category == "Uncharacterized")

list.categories <- as.character(SEED[!duplicated(SEED$Category),1])

SIG.DATAFRAME <- up5
stats_df <- data.frame()
for (i in list.categories) {
  #make dataframe for 1 pathway at a time, consists of pathway (rownames) and gene components
  genegroup <- as.data.frame(subset(SEED, Category == paste(i)))
  genegroup <- genegroup[!duplicated(genegroup$Features),]
  #name of genegroup
  genegroupname<-as.character(genegroup[1,1])
  #get values for our 2x2 contingency table
    a <-nrow(subset(SIG.DATAFRAME, Category == paste(i))) #a is the number of upreg or downreg genes in that PATHWAY 
    b <- nrow(genegroup)#b is the number of other genes in genome (#genes from pathway table - a)
    c <-nrow(SIG.DATAFRAME)-a #c is the number of upreg or downreg genes NOT in that pathway (#sig genes - a)
    d <-nrow(SEED_nodup) #d is the number of other genes in genome not in the pathway 
  #make contingency 2x2 table 
  matrix.stats <- matrix(data=c(a,c,b,d), nrow=2, ncol=2)
  ft<-fisher.test(matrix.stats,alternative = "greater")
  temp1 <- as.data.frame(cbind(ft$p.value, genegroupname))
  stats_df <- rbind(stats_df, temp1)
}
upreg5fisher<-stats_df
colnames(upreg5fisher) <- c("p-value", "genegroup")
## double check output, if there are any zeroes in contingency table, will be reporetd as sig 
write.table(upreg5fisher,"~/GoogleDrive/Stenotrophomonas/output/tables/fishertest5upreg", 
            quote=FALSE, sep = "\t",row.names = F)
write.table(downreg5fisher,"~/GoogleDrive/Stenotrophomonas/output/tables/fishertest5downreg", 
            quote=FALSE, sep = "\t", row.names=F)

####################################################################
### GAGE / PATHVIEW package to identify up or down reg KEGG PATHWAYS
####################################################################

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
