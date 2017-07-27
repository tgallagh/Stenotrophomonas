### 7/7/17 Steno RNA Seq 
#script for my macbook
#R version 3.3.1
#EdgeR analysis of HTSEQcount output
require(dplyr)

setwd("/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/edgeR")
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages("reshape2")
require("edgeR")
require(reshape2) 
counts <- read.table("/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/HTSeq/steno_htseq.xls", header = TRUE)

groups <- c("pH5", "pH5", "pH7", "pH7", "pH7", "pH9", "pH9")
samples <- cbind(groups, colnames(counts[2:8]))
counts.matrix <- counts[,-1]
rownames(counts.matrix) <- counts[,1]
counts.matrix <- counts.matrix[1:4378,]

#edgeR stores data in a simple list-based data object called a DGEList. 
#This type of object is easy to use because it can be manipulated like any list in R. 
#You can make this in R by specifying the counts and the groups in the function DGEList()
d <- DGEList(counts=counts.matrix,group=groups)


#First get rid of genes which did not occur frequently enough. 
#We can choose this cutoff by saying we must have at least 100 counts per million (calculated with cpm() in R) on any particular gene that we want to keep. 
#In this example, we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.
cpm<-(cpm(d))
#returns counts per mILLION FROM A DGElist
#divides raw counts by library size and then multiplies by 1 mil
total_gene<-apply(d$counts, 2, sum) # total gene counts per sample

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep,]
dim(d)

#Recalculate library sizes of the DGEList object
d$samples$lib.size <- colSums(d$counts)
d$samples
#calcNormFactors" function normalizes for RNA composition
#it finds a set of scaling factors for the library sizes (i.e. number of reads in a sample)
#that minimize the log-fold changes between the samples for most genes 
d.norm<- calcNormFactors(d)
#plot MDS - multi-dimensional scaling plot of the RNA samples in which distances correspond to leading log-fold-changes between each pair of RNA samples. 
tiff("bin_NMDS_plot.tiff")
plotMDS<-plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottom", as.character(unique(d$samples$group)), col=1:3, pch=20)
dev.off()
#create a design matrix using model.matrix function
samples <- data.frame(samples)
fac <- paste(samples$groups, sep=".")
fac <- factor(fac)
design <- model.matrix (~0+fac)
colnames(design) <- levels(fac)

#estimate dispersion
d.norm<- estimateGLMCommonDisp(d.norm, design)
d.norm <- estimateGLMTrendedDisp(d.norm, design)
d.norm <- estimateGLMTagwiseDisp(d.norm, design)
#mean variance plot
#grey = raw variance
#light blue = tagwise dispersion varaince
#dark blue = common dispersion
#black line = poission variance
tiff("Meanvarplot_bins.tiff")
meanVarPlot <- plotMeanVar( d.norm, show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )
dev.off()

#BCV plot
#the BCV plot shows the common trended and genewise dispersions
#as a function of average logCPM
#plotBCV(y)

#fit genewise negative binomial general linear model
fit <- glmFit(d.norm, design)

#the likelihood ratio stats are computed for the comparison of interest
#compare (pH 9 group - pH 7 group) to 0:
##### positive logFC for genes upreg in pH 9 compared to pH 7
lrt.7vs9 <- glmLRT(fit, contrast=c(0,-1,1))
#double check coefficiants
#topTags(lrt.7vs9)
#subset79<-subset(P30_P36_RPM, names == "fig|6666666.230262.peg.2616")
#subset79<- melt(subset79)
#plotexample79<- ggplot(data = subset79)
#plotexample79 + geom_point(aes(x=variable, y=value))

#compare (pH 5 group - pH 9 group) to 0:
####positive logFC means higher gene expression in pH 5 > pH 9
lrt.5vs9 <- glmLRT(fit, contrast=c(1,0,-1))
#double check coefficiants
#topTags(lrt.5vs9)
#subset59<-subset(P30_P36_RPM, names == "fig|6666666.230262.peg.4208")
#subset59<- melt(subset59)
#plotexample59<- ggplot(data = subset59)
#plotexample59 + geom_point(aes(x=variable, y=value))

##compare (pH 5 group - pH 7 group) to 0:
####positive logFC means higher gene expression in pH5 
lrt.5vs7 <- glmLRT(fit, contrast=c(1,-1,0))
#######double check the coefficiants by making scatter plots of sig different genes
#topTags(lrt.5vs7)
#subset57<-subset(P30_P36_RPM, names == "fig|6666666.230262.peg.4308")
#subset57 <- melt(subset57)
#plotexample57<- ggplot(data = subset57)
#plotexample57 + geom_point(aes(x=variable, y=value))

#make files with ALL genes
all.5vs7 <- topTags(lrt.5vs7, adjust="BH", n=Inf)
all.5vs9 <- topTags(lrt.5vs9, adjust="BH", n=Inf)
all.7vs9 <- topTags(lrt.5vs7, adjust="BH", n=Inf)
#write out results
write.table(x = all.5vs7, file = "all.5vs7.genes.txt", sep ="\t", quote = FALSE)
#write.table(x = all.5vs9, file = "all.5vs9.BINS.txt", sep ="\t", quote = FALSE)
write.table(x = all.7vs9, file = "all.7vs9.genes.txt", sep ="\t", quote = FALSE)

####FDR < 0.05 (5 in 100 sig genes wrong) AS CUT-OFF
### Write out stats to .csv and also filter by FDR (0.05)
cutoff=0.05

#pH 5 vs pH #7
de.5vs7 <- decideTestsDGE(lrt.5vs7, p=cutoff, adjust="BH")
#obtain number of upreg, downreg genes
summary_5_7<-summary(de.5vs7)
de57tags<-rownames(d)[as.logical(de.5vs7)]
tiff("smearplot_5v7.genes.tiff")
plotSmear(lrt.5vs7, de.tags = de57tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in pH 5 compared to pH 7", cex.main=.8)
dev.off()

#filter by FDR <= cutoff, and logFC (logFC > 1, means 2fold difference)
out5_7 <- topTags(lrt.5vs7, n=Inf, adjust.method="BH")
keep5_7 <- out5_7$table$FDR <= 0.05 & abs(out5_7$table$logFC) >= 1
de5_7<-out5_7[keep5_7,]
write.table(x = de5_7, file = "de5_7_filter.genes.txt", sep ="\t", quote = FALSE)

#pH 7 vs pH 9
de.7vs9 <- decideTestsDGE(lrt.7vs9, p=cutoff, adjust="BH")
summary_7_9<-summary(de.7vs9)
de79tags<-rownames(d)[as.logical(de.7vs9)]
tiff("smearplot_7v9.genes.tiff")
plotSmear(lrt.7vs9, de.tags = de79tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in pH 9 compared to pH 7", cex.main=.8)
dev.off()

out7_9 <- topTags(lrt.7vs9, n=Inf, adjust.method="BH")
keep7_9 <- out7_9$table$FDR <= 0.05 & abs(out7_9$table$logFC) >= 1
de7_9<-out7_9[keep7_9,]
write.table(x = de7_9, file = "de7_9_filter.genes.txt", sep ="\t", quote = FALSE)

#pH 5 vs pH 9
de.5vs9 <- decideTestsDGE(lrt.5vs9, p=cutoff, adjust="BH")
summary_5_9<-summary(de.5vs9)
de59tags<-rownames(d)[as.logical(de.5vs9)]

tiff("smearplot_5v9.BINS.tiff")
plotSmear(lrt.5vs9, de.tags = de59tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in pH 5 compared to pH 9", cex.main=.8)
dev.off()


out5_9 <- topTags(lrt.5vs9, n=Inf, adjust.method="BH")
keep5_9 <- out5_9$table$FDR <= 0.1 & abs(out5_9$table$logFC) >= 1
de5_9<-out5_9[keep5_9,]
write.table(x = de5_9, file = "de5_9_filter.BINS.txt", sep ="\t", quote = FALSE)

#write out summaries from comparisons
summary_all <- cbind(summary_5_7, summary_5_9, summary_7_9)
colnames(summary_all) <- c("5_vs_7", "5_v_9", "7_v_9")
write.table(summary_all, "summary_DGE_all.BINS.txt", sep="\t")

#####################################################
#### edgeR analysis of 1000 bp bins from Steno genome 
#####################################################

alignedreads<-c(3249685, 2411591, 3471746, 2187908, 3607402, 5079716, 4023699)
# ^ these are total # of aligned reads (incldues SE and PE)
working.directory="/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/bowtie2/coverage/"
##############################
### Read in all the data files
##############################
#Coverage by bp from bedtools
setwd(working.directory)
sample.names <- c("P30", "P31", "P32", "P33", "P34", "P35", "P36")
labels <- c("pH5", "pH5", "pH7", "pH7", "pH7", "pH9", "pH9")
sample.info <- cbind(labels, sample.names, alignedreads)
for (file in sample.names) {
  temp <- read.delim(paste(file,"_final.coverage.txt",sep=""),header=FALSE)
  colnames(temp) <- c("scaffold", "position", "coverage")
  assign(file,temp)
}

#GTF file of annotated genes from RAST
GTF <- read.delim(file = "~/GoogleDrive/Stenotrophomonas/data/processed/6666666.230262.gtf",
                  header=FALSE)
GTF$gene <- GTF$V9
GTF <- GTF[-1,]
# rename columns
colnames(GTF) <- c("scaffold", "FIG", "type", "start", "stop", "score",
                   "strand", "frame", "attribute", "gene")
# add a column with just genes
GTF$gene<-gsub(pattern=".*;Name=",replacement= "", x=GTF$gene)




####################################
### Normalize bins by total # counts
####################################
# Calculate counts per million
# CPM = # counts * 10^6 / total # reads 
for (sample in sample.names) {
  reads <- subset(sample.info, sample.names == sample)
  reads <- as.numeric(reads[,3])
  temp <- get(sample)
  temp<- temp %>%
    mutate(CPM =  coverage * 10^6/reads)
  assign(sample,temp)
}
#function for binning:
slidingwindow <- function(windowsize, inputseq){
  starts <- seq(1, length(inputseq-1000), by = 1000)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for(i in 1:n){chunk <- inputseq[starts[i]:(starts[i]+1000-1)]
  chunkmean <- mean(chunk)
  chunkstdv<-sd(chunk)
  chunkbps[i] <- chunkmean
  chunkstats[i]<-chunkstdv}
  return(list(starts,chunkbps,chunkstats))}

### bins for cpm
### bins based off of CPM normalization 
for (sample in sample.names) {
  temp <- get(sample)
  temp$bp <- rownames(temp)
  temp$sample <- sample
  temp.list<-slidingwindow(100, temp[,4])
  df <- data.frame(temp.list[[1]], temp.list[[2]], temp.list[[3]])
  colname<-c("bp","binned_cpm","sd")
  colnames(df)<-colname
  df2 <- merge(x=df, y=temp, by="bp")
  df.name <- paste(sample,"bins", sep="_")
  assign(df.name,df2)
}

ALL.BINS<-cbind(P30_bins$bp, P30_bins$binned_cpm, P31_bins$binned_cpm,
                P32_bins$binned_cpm, P33_bins$binned_cpm, 
                P34_bins$binned_cpm, P35_bins$binned_cpm,
                P36_bins$binned_cpm)
rownames(ALL.BINS) <- ALL.BINS[,1]
ALL.BINS <- ALL.BINS[,-1]
colnames(ALL.BINS) <- c("P30", "P31", "P32",
                        "P33", "P34", "P35", "P36")
ALL.BINS <- na.omit(ALL.BINS)
d.bins <- DGEList(counts=ALL.BINS,group=groups, remove.zeros = TRUE)


#First get rid of genes which did not occur frequently enough. 
#We can choose this cutoff by saying we must have at least 100 counts per million (calculated with cpm() in R) on any particular gene that we want to keep. 
#In this example, we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.
cpm.bins<-(cpm(d.bins))
#returns counts per mILLION FROM A DGElist
#divides raw counts by library size and then multiplies by 1 mil
total_gene<-apply(d.bins$counts, 2, sum) # total gene counts per sample

keep <- rowSums(cpm(d.bins)>1) >= 2
d.bins<- d.bins[keep,]
dim(d.bins)

#Recalculate library sizes of the DGEList object
d.bins$samples$lib.size <- colSums(d.bins$counts)
d.bins$samples
#plot MDS - multi-dimensional scaling plot of the RNA samples in which distances correspond to leading log-fold-changes between each pair of RNA samples. 
setwd("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/")
tiff("bins_NMDS_plot.tiff")
plotMDS<-plotMDS(d.bins, method="bcv", col=as.numeric(d.bins$samples$group))
legend("bottom", as.character(unique(d.bins$samples$group)), col=1:3, pch=20)
dev.off()
#estimate dispersion
d.bins<- estimateGLMCommonDisp(d.bins, design)
d.bins <- estimateGLMTrendedDisp(d.bins, design)
d.bins <- estimateGLMTagwiseDisp(d.bins, design)

#BCV plot
#the BCV plot shows the common trended and genewise dispersions
#as a function of average logCPM
#plotBCV(y)

#fit genewise negative binomial general linear model
fit <- glmFit(d.bins, design)

#the likelihood ratio stats are computed for the comparison of interest
#compare (pH 9 group - pH 7 group) to 0:
##### positive logFC for genes upreg in pH 9 compared to pH 7
lrt.7vs9.bins <- glmLRT(fit, contrast=c(0,-1,1))
#double check coefficiants
#topTags(lrt.7vs9)
#subset79<-subset(P30_P36_RPM, names == "fig|6666666.230262.peg.2616")
#subset79<- melt(subset79)
#plotexample79<- ggplot(data = subset79)
#plotexample79 + geom_point(aes(x=variable, y=value))

#compare (pH 5 group - pH 9 group) to 0:
####positive logFC means higher gene expression in pH 5 > pH 9
lrt.5vs9.bins <- glmLRT(fit, contrast=c(1,0,-1))
#double check coefficiants
#topTags(lrt.5vs9)
#subset59<-subset(P30_P36_RPM, names == "fig|6666666.230262.peg.4208")
#subset59<- melt(subset59)
#plotexample59<- ggplot(data = subset59)
#plotexample59 + geom_point(aes(x=variable, y=value))

##compare (pH 5 group - pH 7 group) to 0:
####positive logFC means higher gene expression in pH5 
lrt.5vs7.bins <- glmLRT(fit, contrast=c(1,-1,0))
#######double check the coefficiants by making scatter plots of sig different genes
#topTags(lrt.5vs7)
#subset57<-subset(P30_P36_RPM, names == "fig|6666666.230262.peg.4308")
#subset57 <- melt(subset57)
#plotexample57<- ggplot(data = subset57)
#plotexample57 + geom_point(aes(x=variable, y=value))

#make files with ALL genes
all.5vs7.bins <- topTags(lrt.5vs7.bins, adjust="BH", n=Inf)
all.5vs9.bins <- topTags(lrt.5vs9.bins, adjust="BH", n=Inf)
all.7vs9.bins <- topTags(lrt.5vs7.bins, adjust="BH", n=Inf)
#write out results
write.table(x = all.5vs7.bins, file = "all.5vs7.bins.txt", sep ="\t", quote = FALSE)
#write.table(x = all.5vs9, file = "all.5vs9.BINS.txt", sep ="\t", quote = FALSE)
write.table(x = all.7vs9.bins, file = "all.7vs9.bins.txt", sep ="\t", quote = FALSE)

####FDR < 0.05 (5 in 100 sig genes wrong) AS CUT-OFF
### Write out stats to .csv and also filter by FDR (0.05)
cutoff=0.05

#pH 5 vs pH #7
de.5vs7.bins <- decideTestsDGE(lrt.5vs7.bins, p=cutoff, adjust="BH")
#obtain number of upreg, downreg genes
summary_5_7<-summary(de.5vs7.bins)
de57tags.bins<-rownames(d.bins)[as.logical(de.5vs7.bins)]
tiff("smearplot_5v7.bins.tiff")
plotSmear(lrt.5vs7.bins, de.tags = de57tags.bins, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated bins in pH 5 compared to pH 7", cex.main=.8)
dev.off()

#filter by FDR <= cutoff, and logFC (logFC > 1, means 2fold difference)
out5_7 <- topTags(lrt.5vs7.bins, n=Inf, adjust.method="BH")
keep5_7 <- out5_7$table$FDR <= 0.05 & abs(out5_7$table$logFC) >= 1
de5_7.bins<-out5_7[keep5_7,]
write.table(x = de5_7.bins, file = "de5_7_filter.bins.txt", sep ="\t", quote = FALSE)

#pH 7 vs pH 9
de.7vs9.bins <- decideTestsDGE(lrt.7vs9.bins, p=cutoff, adjust="BH")
summary_7_9<-summary(de.7vs9)
de79tags<-rownames(d)[as.logical(de.7vs9)]
tiff("smearplot_7v9.bins.tiff")
plotSmear(lrt.7vs9.bins, de.tags = de79tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated bins in pH 9 compared to pH 7", cex.main=.8)
dev.off()

out7_9 <- topTags(lrt.7vs9.bins, n=Inf, adjust.method="BH")
keep7_9 <- out7_9$table$FDR <= 0.05 & abs(out7_9$table$logFC) >= 1
de7_9.bins<-out7_9[keep7_9,]
write.table(x = de7_9.bins, file = "de7_9_filter.bins.txt", sep ="\t", quote = FALSE)

de7_9.bins.df<-data.frame(de7_9.bins)
de5_7.bins.df<-data.frame(de5_7.bins)

BP<-as.numeric(rownames(de5_7.bins.df))
BP <- BP[2]
SUBSET<-subset(P30_bins, bp== BP)
SCAFFOLD <- as.character(SUBSET$scaffold)
SCAFFOLD.POSITION.START <- as.numeric(SUBSET$position)
SCAFFOLD.POSITION.STOP <- SCAFFOLD.POSITION.START + 999
SUBSET.GTF <- GTF %>%
  filter(scaffold == SCAFFOLD) %>%
  filter(start > SCAFFOLD.POSITION.START |
           start > SCAFFOLD.POSITION.STOP) %>%
  filter(stop < SCAFFOLD.POSITION.STOP)
         
         |
           stop > SCAFFOLD.POSITION.START)
  
  | 
          start < SCAFFOLD.POSITION.STOP)


