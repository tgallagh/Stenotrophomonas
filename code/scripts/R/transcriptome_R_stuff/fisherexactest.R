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

up9 <- subset(data,which=="9up")
up9 <- subset(up9, !Category == "Uncharacterized")

down9 <- subset(data,which=="9down")
down9 <- subset(down9, !Category == "Uncharacterized")


list.categories <- as.character(SEED[!duplicated(SEED$Subcategory),1])
list.categories <- as.character(unique(SEED$Subcategory))

SIG.DATAFRAME <- down5
stats_df <- data.frame()
for (i in list.categories) {
  #make dataframe for 1 pathway at a time, consists of pathway (rownames) and gene components
  genegroup <- as.data.frame(subset(SEED, Subcategory == paste(i)))
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

