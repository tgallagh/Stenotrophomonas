#### Steno RNA seq
#### updated 05/18/17
#### R version 3.3.1
require(dplyr)
require(ggplot2)
require(reshape2)

### Assign labels and info to variables 
sample.names <- c("P30", "P31", "P32", "P33", "P34", "P35", "P36")
labels <- c("pH5", "pH5", "pH7", "pH7", "pH7", "pH9", "pH9")
alignedreads<-c(3249685, 2411591, 3471746, 2187908, 3607402, 5079716, 4023699)
# ^ these are total # of aligned reads (incldues SE and PE)
sample.info <- cbind(labels, sample.names, alignedreads)
working.directory="/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/bowtie2/coverage/"
##############################
### Read in all the data files
##############################
#Coverage by bp from bedtools
setwd(working.directory)
for (file in sample.names) {
  temp <- read.delim(paste(file,"_final.coverage.txt",sep=""),header=FALSE)
  colnames(temp) <- c("scaffold", "position", "coverage")
  assign(file,temp)
}

#GTF file of annotated genes from RAST
GTF <- read.delim(file = "~/GoogleDrive/Stenotrophomonas/data/processed/6666666.230262.gtf",
                  header=FALSE)
GTF$gene <- GTF$V9
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

##############################
# Graphs of genes of interest 
#############################
all_bp_cov <- cbind(as.data.frame(P30$scaffold), P30$position, P30$coverage,
                    P31$coverage, P32$coverage, P33$coverage, P34$coverage,
                    P35$coverage, P36$coverage)

# get number of basepairs for whole genome 
all_bp_cov$bp <- 1:nrow(all_bp_cov)
#rename columns 
colnames(all_bp_cov) <- c("scaffold", "position", "P30", 
                          "P31" , "P32", "P33", "P34", "P35" , "P36", "bp")

ALL.COV = all_bp_cov

#make color panel to color lines for plots below 
sample.colors<- c("P30"="red", "P31" = "red" ,
                  "P32" = "black", "P33" = "black",
                  "P34" = "black", "P35" = "blue" ,
                  "P36" = "blue") 


### this function will give single bp resolution for a gene of interest

make.gene.plots <- function(x) {
GENE.INFO <- subset(GTF, gene == x)
GENE.START <- GENE.INFO$start
GENE.STOP <- GENE.INFO$stop 
GENE.SCAFFOLD <- as.character(GENE.INFO$scaffold)
SAMPLE.GENE <- subset(ALL.COV, scaffold ==GENE.SCAFFOLD)
SAMPLE.GENE <- SAMPLE.GENE %>%
  filter(position > GENE.START & position < GENE.STOP)
plot <- ggplot()
plot + geom_line(data=SAMPLE.GENE, aes(x=bp, y=P30, color="P30"))+
  geom_line(data=SAMPLE.GENE, aes(x=bp, y=P31, color="P31"))+
  geom_line(data=SAMPLE.GENE, aes(x=bp, y=P32, color="P32"))+
  geom_line(data=SAMPLE.GENE, aes(x=bp, y=P33, color="P33"))+
  geom_line(data=SAMPLE.GENE, aes(x=bp, y=P34, color="P34"))+
  geom_line(data=SAMPLE.GENE, aes(x=bp, y=P35, color="P35"))+
  geom_line(data=SAMPLE.GENE, aes(x=bp, y=P36, color="P36"))+
                   scale_color_manual(values=sample.colors, name="Sample")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Counts per million") +
  xlab("Genome position")+
  ggtitle(paste(x))
}

### make.gene.plots(x) <--- change x to put gene of interest 
x = "DNA gyrase subunit A (EC 5.99.1.3);Ontology_term=KEGG_ENZYME:5.99.1.3"
make.gene.plots(x)
x = GTF[2727,10] #"reca"
make.gene.plots(x)
x = GTF[3924,10] #transcription termination factor rho
make.gene.plots(x)
x = GTF[3572,10] #"DNA directed RNA polymerase Beta subunit"
make.gene.plots(x) 
x = GTF[3864,10] #"DNA directed RNA polymerase Beta subunit"
make.gene.plots(x)

#pick housekeeping genes based on this criteria:
  # CPM across the gene > 1
  # list of housekeeping genes: https://www.ncbi.nlm.nih.gov/pubmed/26149127
# determine change in average coverage for those housekeeping genes

calculate.coverage <- function(HKG) {
  GENE.INFO <- subset(GTF, gene == HKG)
  GENE.START <- GENE.INFO$start
  GENE.STOP <- GENE.INFO$stop 
  GENE.SCAFFOLD <- as.character(GENE.INFO$scaffold)
  SAMPLE.GENE <- subset(ALL.COV, scaffold ==GENE.SCAFFOLD)
  SAMPLE.GENE <- SAMPLE.GENE %>%
    filter(position > GENE.START & position < GENE.STOP)
  PH5.AVERAGE <- mean((SAMPLE.GENE$P30 + SAMPLE.GENE$P31)/2)
  PH7.AVERAGE <- mean((SAMPLE.GENE$P32 + SAMPLE.GENE$P33 +
                         SAMPLE.GENE$P34)/3)
  PH9.AVERAGE <- mean((SAMPLE.GENE$P35 + SAMPLE.GENE$P36)/2)
  HKG.AVERAGE <- cbind(PH5.AVERAGE, PH7.AVERAGE, PH9.AVERAGE)
  assign("HKG.avg",HKG.AVERAGE, envir = .GlobalEnv)
}
HKG = "DNA gyrase subunit A (EC 5.99.1.3);Ontology_term=KEGG_ENZYME:5.99.1.3"
calculate.coverage(HKG)
housekeeping <- HKG.avg
HKG = GTF[2727,10] #"reca"
calculate.coverage(HKG)
housekeeping <- rbind(housekeeping, HKG.avg)
HKG = GTF[3924,10] #transcription termination factor rho
calculate.coverage(HKG)
housekeeping <- rbind(housekeeping, HKG.avg)
HKG = GTF[3572,10] #"DNA directed RNA polymerase Beta subunit"
calculate.coverage(HKG)
housekeeping <- rbind(housekeeping, HKG.avg)
HKG = GTF[3864,10] #"DNA directed RNA polymerase Beta subunit"
calculate.coverage(HKG)
housekeeping <- rbind(housekeeping, HKG.avg)
housekeeping<-colMeans(housekeeping)/100
pH9.change <- housekeeping[3] / housekeeping[2]
pH5.change <- housekeeping[1] / housekeeping[2]
pH7.change <- 1 
housekeeping.final<- cbind(pH5.change, pH7.change, pH9.change)

###############################
################################
#Make 100 bp bins 
################################
#################################

#function for binning:
slidingwindow <- function(windowsize, inputseq){
  starts <- seq(1, length(inputseq-100), by = 100)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for(i in 1:n){chunk <- inputseq[starts[i]:(starts[i]+100-1)]
  chunkmean <- mean(chunk)
  chunkstdv<-sd(chunk)
  chunkbps[i] <- chunkmean
  chunkstats[i]<-chunkstdv}
  return(list(starts,chunkbps,chunkstats))}

#loop over samples

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


### bins for coverage 
### bins based off of non-normalized coverage 
#### this should have non-normalized bp coverage in 3rd column of dataframe

for (sample in sample.names) {
  temp <- get(sample)
  temp$bp <- rownames(temp)
  temp$sample <- sample
  temp.list<-slidingwindow(100, temp[,3])
  df <- data.frame(temp.list[[1]], temp.list[[2]], temp.list[[3]])
  colname<-c("bp","binned_cov","sd")
  colnames(df)<-colname
  df2 <- merge(x=df, y=temp, by="bp")
  df.name <- paste(sample,"bins_cov", sep="_")
  assign(df.name,df2)
}

# make new dataframes that take into account housekeeping gene normalization
P30_bins<- P30_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[1])
P31_bins<- P31_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[1])
P32_bins<- P32_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[2])
P33_bins<- P33_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[2])
P34_bins<- P34_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[2])
P35_bins<- P35_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[3])
P36_bins<- P36_bins %>%
  mutate(hkg_norm = binned_cpm / housekeeping.final[3])


#### calculate average coverage for each sample
#### put dataframes into list to use apply function 
LIST <- list(P30_bins, P31_bins, P32_bins,
               P33_bins, P34_bins, P35_bins, 
               P36_bins)
#### Function for calculating average CPM across genome 
average.coverage <- function(INPUT,OUTPUT){
  cpm.col <- as.numeric(INPUT$binned_cpm)
  cpm.col.mean <- mean(cpm.col, na.rm=TRUE)
  assign(paste(OUTPUT),cpm.col.mean, envir = .GlobalEnv)
}

#apply to list of dataframes
allcoverageaverage<-lapply(LIST, FUN= average.coverage, OUTPUT="cov")
allcoverageaverage <- unlist(allcoverageaverage)
allcoverageaverage <- cbind(allcoverageaverage, sample.names)

#### Function for calculating average coverage using non-normalized bp coverage
average.coverage.raw <- function(INPUT,OUTPUT){
  cpm.col <- as.numeric(INPUT$binned_cov)
  cpm.col.mean <- mean(cpm.col, na.rm=TRUE)
  assign(paste(OUTPUT),cpm.col.mean, envir = .GlobalEnv)
}

LIST2 <- list(P30_bins_cov, P31_bins_cov, P32_bins_cov,
             P33_bins_cov, P34_bins_cov, P35_bins_cov, 
             P36_bins_cov)

allcoverageaverage.bp<-lapply(LIST2, FUN= average.coverage.raw, OUTPUT="cov")
allcoverageaverage.bp <- unlist(allcoverageaverage.bp)
allcoverageaverage.bp <- cbind(allcoverageaverage.bp, sample.names)


### normalize CPM bins by average coverage of CPM normalized 
#### this is weird, don't know why I did it 

P30_bins<- P30_bins %>%
  mutate(avg_norm = binned_cpm * 30 / as.numeric(allcoverageaverage[1,1]))

P31_bins<- P31_bins %>%
  mutate(avg_norm = binned_cpm * 30 / as.numeric(allcoverageaverage[2,1]))

P32_bins<- P32_bins %>%
  mutate(avg_norm = binned_cpm * 30 / as.numeric(allcoverageaverage[3,1]))

P33_bins<- P33_bins %>%
  mutate(avg_norm = binned_cpm * 30 / as.numeric(allcoverageaverage[4,1]))

P34_bins<- P34_bins %>%
  mutate(avg_norm = binned_cpm * 30/ as.numeric(allcoverageaverage[5,1]))

P35_bins<- P35_bins %>%
  mutate(avg_norm = binned_cpm * 30 / as.numeric(allcoverageaverage[6,1]))

P36_bins<- P36_bins %>%
  mutate(avg_norm = binned_cpm * 30 / as.numeric(allcoverageaverage[7,1]))


#### calculate mean basepair coverage of all samples
mean(as.numeric(allcoverageaverage.bp[,1]))

P30_bins_cov<- P30_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[1,1]))

P31_bins_cov<- P31_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[2,1]))

P32_bins_cov<- P32_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[3,1]))

P33_bins_cov<- P33_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[4,1]))

P34_bins_cov<- P34_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[5,1]))

P35_bins_cov<- P35_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[6,1]))

P36_bins_cov<- P36_bins_cov %>%
  mutate(avg_norm = binned_cov * 119 / as.numeric(allcoverageaverage.bp[7,1]))

#### plot individual samples
plot <- ggplot()
#color panel for different conditions 
cols <- c("ph5"="red","ph7"="black","ph9"="blue")
YVALUE = "log2(avg_norm)"


######## GLOBAL PLOT FOR RNASEQ 
plot + geom_line(data=P30_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph5")), size=4)+
  scale_color_manual(values=cols, name="Condition")+
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14, colour="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("counts per million")+
  xlab("Genome Position") +
  ggtitle("Coverage of all samples normalized by \n mean genome-wide coverage")+
  geom_line(data=P31_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph5")), size=4)+
   geom_line(data=P32_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph7")),  size=3)+
  geom_line(data=P33_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph7")), size=3)+
  geom_line(data=P34_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph7")), size=3)+
  geom_line(data=P35_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph9")), size=1)+
  geom_line(data=P36_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph9")), size=1)
  
ylim(0,300000)

cols2 <-  c("ph5"="pink","ph7"="gray","ph9"="#42cbf4", 
            "ph5_hkg" = "red", "ph7_hkg" = "black", 
            "ph9_hkg" = "blue")

plot + 
  geom_line(data=P30_bins, aes(x=bp, y=binned_cpm, color="ph5"), size=5)+
  geom_line(data=P31_bins, aes(x=bp, y=binned_cpm,color="ph5"), size=5)+
  geom_line(data=P32_bins, aes(x=bp, y=binned_cpm,color="ph7"), size=4)+
  geom_line(data=P33_bins, aes(x=bp, y=binned_cpm,color="ph7"), size=4)+
  geom_line(data=P34_bins, aes(x=bp, y=binned_cpm,color="ph7"), size=4)+
  #geom_line(data=P35_bins, aes(x=bp, y=binned_cpm,color="ph9"), size=3)+
  #geom_line(data=P36_bins, aes(x=bp, y=binned_cpm,color="ph9"), size=3)+
  geom_line(data=P30_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph5_hkg")), size=2)+
  scale_color_manual(values=cols2, name="Condition")+
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14, colour="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("counts per million")+
  xlab("Genome Position") +
  ggtitle("Comparison of normalization steps")+
  geom_line(data=P31_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph5_hkg")), size=2)+
  geom_line(data=P32_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph7_hkg")),  size=1.5)+
  geom_line(data=P33_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph7_hkg")), size=1.5)+
  geom_line(data=P34_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph7_hkg")), size=1.5)+
 # geom_line(data=P35_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph9_hkg")), size=1)+
  #geom_line(data=P36_bins, aes_string(x="bp", y=YVALUE, color=shQuote("ph9_hkg")), size=1)+
  ylim(0,300000)
  

#### plot of non-normalized binned bp coverage

YVALUE = "avg_norm"
plot + geom_line(data=P30_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph5")), size=4)+
  scale_color_manual(values=cols, name="Condition")+
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14, colour="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("counts per million")+
  xlab("Genome Position") +
  ggtitle("COVERAGE OF RNA READS ACROSS GENOME
          NORMALIZED BY GENOME-WIDE EXPRESSION")+
  geom_line(data=P31_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph5")), size=4)+
  geom_line(data=P32_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph7")),  size=3)+
  geom_line(data=P33_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph7")), size=3)+
  geom_line(data=P34_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph7")), size=3)+
  geom_line(data=P35_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph9")), size=1)+
  geom_line(data=P36_bins_cov, aes_string(x="bp", y=YVALUE, color=shQuote("ph9")), size=1)+
  ylim(0,300000)

#log2 CPM
plot + 
  geom_line(data=P30_bins, aes(x=bp, y=log2(binned_cpm), color="ph5"), size=1.5)+
  geom_line(data=P31_bins, aes(x=bp, y=log2(binned_cpm), color="ph5"), size=1.5)+
  geom_line(data=P32_bins, aes(x=bp, y=log2(binned_cpm), color="ph7"), size=1)+
  geom_line(data=P33_bins, aes(x=bp, y=log2(binned_cpm), color="ph7"), size=1)+
  geom_line(data=P34_bins, aes(x=bp, y=log2(binned_cpm), color="ph7"), size=1)+
  geom_line(data=P35_bins, aes(x=bp, y=log2(binned_cpm), color="ph9"), size=1)+
  geom_line(data=P36_bins, aes(x=bp, y=log2(binned_cpm), color="ph9"), size=1)+
  scale_color_manual(values=cols, name="Condition")+
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14, colour="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("Log10 Counts per million")+
  xlab("Genome Position")


#### put all binning data into one data frame 
all_bins <- rbind(P30_bins, P31_bins, P32_bins, P33_bins,P34_bins, P35_bins, P36_bins)
all_bins <- all_bins[,-(6:7)]
all_bins.m <- melt(all_bins , id.vars=c("bp", "scaffold", "sample", "position"),
                   measure.vars=c("binned_cpm", "sd"))

## calculate averages for 3 conditions 
pH5_avg <- as.data.frame(P30_bins$binned_cpm+P31_bins$binned_cpm)/2
pH7_avg <- as.data.frame(P32_bins$binned_cpm+P33_bins$binned_cpm + P34_bins$binned_cpm)/3
pH9_avg <- as.data.frame(P35_bins$binned_cpm + P36_bins$binned_cpm)/2
averages_bins <- cbind(pH5_avg, pH7_avg, pH9_avg)
averages_bins <- cbind(averages_bins, P30_bins$bp, P30_bins$scaffold, P30_bins$position)
colnames(averages_bins)<- c("pH5", "pH7","pH9", "bp", "scaffold", "scaffold.position")

pH5_avg.hkg <- as.data.frame(P30_bins$hkg_norm+P31_bins$hkg_norm)/2
pH7_avg.hkg <- as.data.frame(P32_bins$hkg_norm+P33_bins$hkg_norm + P34_bins$hkg_norm)/3
pH9_avg.hkg <- as.data.frame(P35_bins$hkg_norm + P36_bins$hkg_norm)/2
averages_bins <- cbind(pH5_avg, pH7_avg, pH9_avg)
averages_bins <- cbind(averages_bins, P30_bins$bp, P30_bins$scaffold, P30_bins$position)
averages_bins <- cbind(averages_bins, pH5_avg.hkg, pH7_avg.hkg, pH9_avg.hkg  )
colnames(averages_bins)<- c("pH5", "pH7","pH9", "bp", "scaffold", "scaffold.position", 
                            "pH5.hkg", "pH7.hkg", "pH9.hkg")


pH5_avg.gw <- as.data.frame(P30_bins_cov$avg_norm+P31_bins_cov$avg_norm)/2
pH7_avg.gw <- as.data.frame(P32_bins_cov$avg_norm+P33_bins_cov$avg_norm + P34_bins_cov$avg_norm)/3
pH9_avg.gw<- as.data.frame(P35_bins_cov$avg_norm + P36_bins_cov$avg_norm)/2
averages_bins.gw <- cbind(pH5_avg.gw, pH7_avg.gw, pH9_avg.gw)
averages_bins.gw <- cbind(averages_bins.gw, P30_bins$bp, P30_bins$scaffold, P30_bins$position)

colnames(averages_bins.gw)<- c("pH5.gw", "pH7.gw","pH9.gw", "bp", "scaffold", "scaffold.position")

## calculate log2FC for different groups
averages_bins <- averages_bins %>%
  mutate(log2fc_ph57 = log2(pH5/pH7) )
averages_bins <- averages_bins %>%
  mutate(log2fc_ph97 = log2(pH9/pH7) )
averages_bins <- averages_bins %>%
  mutate(log2fc_ph59 = log2(pH5/pH9) )

averages_bins <- averages_bins %>%
  mutate(log2fc_ph57.hkg = log2(pH5.hkg/pH7.hkg) )
averages_bins <- averages_bins %>%
  mutate(log2fc_ph97.hkg = log2(pH9.hkg/pH7.hkg) )
averages_bins <- averages_bins %>%
  mutate(log2fc_ph59.hkg = log2(pH5.hkg/pH9.hkg) )


averages_bins.gw <- averages_bins.gw %>%
  mutate(log2fc_ph57.gw = log2(pH5.gw/pH7.gw) )
averages_bins.gw <- averages_bins.gw %>%
  mutate(log2fc_ph97.gw = log2(pH9.gw/pH7.gw) )
averages_bins.gw <- averages_bins.gw %>%
  mutate(log2fc_ph59.gw = log2(pH5.gw/pH9.gw) )

averages_bins_all <- cbind(averages_bins, averages_bins.gw)
averages_bins_all <- averages_bins_all[,-(13:15)]

## Filter by logFC 
## log2(fc) > 1 or < -1
## filter out low counts (only bins with > 1)


pH5_pH7 <- averages_bins_all %>%
  filter(pH5 > 1) %>% 
  filter(pH7 >1) %>%
  filter(log2fc_ph57.gw > 1 | log2fc_ph57.gw < -1) 

#pH5_pH7.hkg <- averages_bins %>%
  filter(pH5 > 1) %>% 
  filter(pH7 > 1) %>%
  filter(log2fc_ph57.hkg > 1 | log2fc_ph57.hkg < -1) 

pH7_pH9 <- averages_bins %>%
  filter(pH7 > 1) %>% 
  filter(pH9 > 1) %>%
  filter(log2fc_ph97 > 1 | log2fc_ph97 < -1) 

#pH7_pH9.hkg <- averages_bins %>%
  filter(pH7 > 1) %>%   #### filter by cpm 
  filter(pH9 > 1) %>%
  filter(log2fc_ph97.hkg > 1 | log2fc_ph97.hkg < -1) 


################################
######## Operon identification
###############################

## function to identify operons (maybe?)
## set number of bins to 20 (2000 bp up or down regulated)

### this is the function
### input dataframe must have bp in 4th column


sig.operons.function <- function(df, number_bins, df_name) {
#set variables 
  PREVIOUS.POSITION = 0
  COUNTER = 1
  DF_OUTPUT <- c()
  TEMP.OPERON <- c()
  n=0
# ENTER LOOP 
    for (i in 1:nrow(df)) {
      # loop through rows that represent bin number in coverage data frame
      # If the suceeding row is the next bin or the 2nd bin 
           # then add 1 to counter
      # allows for 100 bp bins that do not follow trend of up or down reg
      if (df[i,4] == (PREVIOUS.POSITION + 100) 
          | df[i,4] == (PREVIOUS.POSITION + 200)) {
        COUNTER = COUNTER + 1
      }  else { # exit the counter once we fail to get succeeding bins
        if (COUNTER > number_bins) {
          #if number of bins exceeds user threshold 
          # then make data frame with info about operons
          # assign bins to operon number 
          n= n+1
          TEMP.OPERON <- df[(i-(COUNTER-1)):i-1,]
          TEMP.OPERON <- TEMP.OPERON %>%
            mutate(operon_no = n)
          df_output <- rbind(df_output, TEMP.OPERON) ## append to dataframe
        }
       COUNTER=1 #reset counter to 1 
      }
      previous.position=df[i,4]
    } # EXIT LOOP 
    assign(paste(df_name),df_output, envir = .GlobalEnv)
    }

### user inputs input dataframe, number of bins (20 = 2000 bp limit), and output df name

sig.operons.function(pH5_pH7, 10, "pH5_pH7_operons.gw")
sig.operons.function(pH7_pH9.hkg, 9, "pH9_pH7_operons.hkg")

###### 
###### filter gtf table so that we get genes falling within this operon

get.annotation.info <- function(df, number) {
  data=subset(df, operon_no == paste(number))
  scaffold.position<-data$scaffold.position
  scaffold.position.first <- scaffold.position[1]
  scaffold.position.last <- scaffold.position[length(scaffold.position)]
  operon.length <- abs(scaffold.position.last - scaffold.position.first)
  SCAFFOLD <- as.character(data$scaffold[1])
  subsetgtf <- subset(GTF, scaffold == SCAFFOLD)
  subsetgtf<- subsetgtf %>% filter (start > (scaffold.position.first - (operon.length/2))
                                    & stop < scaffold.position.last + (operon.length/2))
  # convert scaffold number to genome position to fit onto the graph 
  CONVERT.FACTOR <- data$bp[1] - data$scaffold.position[1]
  subsetgtf <- subsetgtf %>% 
    mutate(genome.start = round(CONVERT.FACTOR + start))
  subsetgtf <- subsetgtf %>% 
    mutate(genome.stop = round(CONVERT.FACTOR + stop))
assign("data", data, envir =  .GlobalEnv)
assign("subsetgtf", subsetgtf, envir=.GlobalEnv)
}

get.annotation.info(pH5_pH7_operons.gw, 4)

  plot <- ggplot()
    plot + geom_line(data=data, aes(x=bp, y=pH5, color="ph5"))+
      geom_line(data=data, aes(x=bp, y=pH7, color="ph7"), size=1) +
      geom_line(data=data, aes(x=bp, y=pH9, color="ph9"), size=1) +
      theme(axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=14, colour="black"),
            axis.title.y=element_text(size=16),
            axis.title.x=element_text(size=14, colour="black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
      scale_color_manual(values=cols, name="Condition") + 
  #geom_segment()+
  xlab("Genome Position") +
  ylab("Normalized Coverage")+
  ggtitle("Downregulated operon in pH5")+
  geom_segment(data=subsetgtf[1,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4)+
      geom_segment(data=subsetgtf[2,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="green",size=4)+
    geom_segment(data=subsetgtf[3,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4)+
    #  geom_segment(data=subsetgtf[4,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4) +
    #  geom_segment(data=subsetgtf[5,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="green",size=4)+
      #geom_segment(data=subsetgtf[6,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4)+
 #geom_label(data=subsetgtf[1,], aes(label=gene, x=(genome.start+1500),y=-5), size=3, fill="purple")+
    geom_label(data=subsetgtf[2,], aes(label=gene, x=(genome.start+1800),y=-5), size=3, fill="green")
  #  geom_label(data=subsetgtf[4,], aes(label=gene, x=genome.start,y=-15), size=3, fill="purple")
    
      
      
    
    
  ### list of upregulated genes:
      #### alkyl hyperoixde reductase protein C
      ###  
     +
  
      
      geom_segment(data=subsetgtf[2,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="green",
               size=4)+
  geom_label(data=subsetgtf[2,], aes(label=gene, x=982000,y=-50), fill="green", size=3) +
  geom_segment(data=subsetgtf[3,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",
               size=4)+
  geom_label(data=subsetgtf[3,], aes(label=gene, x=genome.start,y=-75), fill="purple", size=3) +
  geom_segment(data=subsetgtf[4,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="green",
               size=4)+
  geom_label(data=subsetgtf[4,], aes(label=gene, x=genome.start,y=-100), fill="green", size=3) +
  geom_segment(data=subsetgtf[5,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",
               size=4)+
geom_label(data=subsetgtf[5,], aes(label=gene, x=genome.start,y=-125), fill="purple", size=3) 



plot <- ggplot()
plot + geom_point(data=data, aes(x=bp, y=log2fc_ph57))
                  
                  
  geom_line(data=data, aes(x=bp, y=pH7, color="ph7"), size=1) +
  geom_line(data=data, aes(x=bp, y=pH9, color="ph9"), size=1) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14, colour="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  scale_color_manual(values=cols, name="Condition") + 
  #geom_segment()+
  ylab("Counts per million")+
  xlab("Genome Position")


