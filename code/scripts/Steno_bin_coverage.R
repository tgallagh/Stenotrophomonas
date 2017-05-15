#### Steno RNA seq
#### updated 05/15/17
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
GTF <- read.delim(file = "/Users/Tara/GoogleDrive/P1_RNAseq/data/processed/annotations/6666666.246352.gtf",
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
x = "Alginate biosynthesis protein AlgJ"
make.gene.plots(x)

###############################
################################
#Make 100 bp bins 
################################
#################################

#function for binning:
slidingwindow <- function(windowsize, inputseq){
  starts <- seq(1, length(P30[,3])-100, by = 100)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for(i in 1:n){chunk <- P30[,3][starts[i]:(starts[i]+100-1)]
  chunkmean <- mean(chunk)
  chunkstdv<-sd(chunk)
  chunkbps[i] <- chunkmean
  chunkstats[i]<-chunkstdv}
  return(list(starts,chunkbps,chunkstats))}

#loop over samples
for (sample in sample.names) {
  temp <- get(sample)
  temp$bp <- rownames(temp)
  temp$sample <- sample
  temp.list<-slidingwindow(100, temp[,3])
  df <- data.frame(temp.list[[1]], temp.list[[2]], temp.list[[3]])
  colname<-c("bp","binned_cpm","sd")
  colnames(df)<-colname
  df2 <- merge(x=df, y=temp, by="bp")
  df.name <- paste(sample,"bins", sep="_")
  assign(df.name,df2)
}

#### plot individual samples
plot <- ggplot()
#color panel for different conditions 
cols <- c("ph5"="red","ph7"="black","ph9"="blue")
plot + geom_line(data=P30_bins, aes(x=bp, y=binned_cpm, color="ph5"), size=4)+
  geom_line(data=P31_bins, aes(x=bp, y=binned_cpm, color="ph5"), size=4)+
  geom_line(data=P32_bins, aes(x=bp, y=binned_cpm, color="ph7"), size=3)+
  geom_line(data=P33_bins, aes(x=bp, y=binned_cpm, color="ph7"), size=3)+
  geom_line(data=P34_bins, aes(x=bp, y=binned_cpm, color="ph7"), size=3)+
  geom_line(data=P35_bins, aes(x=bp, y=binned_cpm, color="ph9"), size=1)+
  geom_line(data=P36_bins, aes(x=bp, y=binned_cpm, color="ph9"), size=1)+
  scale_color_manual(values=cols, name="Condition")+
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14, colour="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("Counts per million")+
  xlab("Genome Position")

#log2 CPM
plot + geom_line(data=P30_bins, aes(x=bp, y=log2(binned_cpm), color="ph5"), size=4)+
  geom_line(data=P31_bins, aes(x=bp, y=log2(binned_cpm), color="ph5"), size=4)+
  geom_line(data=P32_bins, aes(x=bp, y=log2(binned_cpm), color="ph7"), size=3)+
  geom_line(data=P33_bins, aes(x=bp, y=log2(binned_cpm), color="ph7"), size=3)+
  geom_line(data=P34_bins, aes(x=bp, y=log2(binned_cpm), color="ph7"), size=3)+
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

## calculate log2FC for different groups
averages_bins <- averages_bins %>%
  mutate(log2fc_ph57 = log2(pH5/pH7) )
averages_bins <- averages_bins %>%
  mutate(log2fc_ph97 = log2(pH9/pH7) )
averages_bins <- averages_bins %>%
  mutate(log2fc_ph59 = log2(pH5/pH9) )

## Filter by logFC 
## log2(fc) > 1 or < -1
## filter out low counts (only bins with > 1)
pH5_pH7 <- averages_bins[,1:7] %>%
  filter(pH5 > 1) %>% 
  filter(pH7 >1) %>%
  filter(log2fc_ph57 > 1 | log2fc_ph57 < -1) 

################################
######## Operon identification
###############################

## function to identify operons (maybe?)
## set number of bins to 20 (2000 bp up or down regulated)

### this is the function
sig.operons.function <- function(df, number_bins, df_name) {
previous.position = 0
counter = 0
df_output <- c()
temp.operon <- c()
n=0
df = pH5_pH7
for (i in 1:nrow(df)) {
  if (df[i,4] == previous.position + 100) {
    counter = counter + 1
  } else {
    if (counter > number_bins) {
      n= n+1
      temp.operon <- df[(i-(counter-1)):i-1,]
      temp.operon <- temp.operon %>%
        mutate(operon_no = n)
      df_output <- rbind(df_output, temp.operon)
    }
    counter=1
  }
  previous.position=df[i,4]
}
assign(paste(df_name),df_output, envir=.GlobalEnv)
}

### user inputs input dataframe, number of bins (20 = 2000 bp limit), and output df name
sig.operons.function(pH5_pH7, 20, "pH5_pH7_operons")

###### 
###### filter gtf table so that we get genes falling within this operon

data=subset(pH5_pH7_operons, operon_no == "1")
scaffold.position<-data$scaffold.position
scaffold.position.1 <- scaffold.position[1]
scaffold.position.last <- scaffold.position[length(scaffold.position)]
subsetgtf <- subset(GTF, scaffold == "scaffold3.1")
subsetgtf<- subsetgtf %>% filter (start > scaffold.position.1 & stop < scaffold.position.last)
# convert scaffold number to genome position to fit onto the graph 
subsetgtf <- subsetgtf %>% 
  mutate(genome.start = round(((980201-67192) + start)))
subsetgtf <- subsetgtf %>% 
  mutate(genome.stop = round(((980201-67192) + stop)))

######
plot <- ggplot()
plot + geom_line(data=data, aes(x=bp, y=pH5, color="ph5"))+
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
  xlab("Genome Position") +
  ylab("Counts per million")+
  ggtitle("Downregulated 'operon' in pH5")+
  geom_segment(data=subsetgtf[1,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",
               size=4)+
  geom_label(data=subsetgtf[1,], aes(label=gene, x=genome.start,y=-25), fill="purple", size=3) +
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


