### Operon Identificaton at single bp resolution for RNAseq data
#### TARA
##### R 3.3.1
###### Updated: 7/19/17


#### Load libraries
require(dplyr)
require(reshape2)
require(xlsx)
require(ggplot2)

# name of individual samples 
SAMPLE.NAMES <- c("P30", "P31", "P32", "P33", "P34", "P35", "P36") 
# conditions of samples in same order as SAMPLE.NAMES vector
SAMPLE.LABELS<- c("pH5", "pH5", "pH7", "pH7", "pH7", "pH9", "pH9")
#Total number of aligned reads, in same order as SAMPLE.NAMES vector
ALIGNED.READS<-c(3249685, 2411591, 3471746, 2187908, 3607402, 5079716, 4023699)
SAMPLE.INFO <- cbind(SAMPLE.NAMES,SAMPLE.LABELS, ALIGNED.READS)
#set working directory with bedtools coverage file 
WORKING.DIRECTORY= "/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/bowtie2/coverage/"

##############################
### Read in data files
##############################
setwd(WORKING.DIRECTORY)

# for loop for reading in all data files with sample name listed in sample.names
# note, suffix of coverage file name must be "_final.coverage.txt" 
# and prefix must be sample name
for (file in SAMPLE.NAMES) {
  temp <- read.delim(paste(file,"_final.coverage.txt",sep=""),header=FALSE)
  colnames(temp) <- c("scaffold", "position", "coverage")
  assign(file,temp)
}
# Read in GTF file of annotated genes from RAST 
GTF <- read.delim(file = "~/GoogleDrive/Stenotrophomonas/data/processed/6666666.230262.gtf",
                  header=FALSE)
GTF$gene <- GTF$V9
#rename columns 
colnames(GTF) <- c("scaffold", "FIG", "type", "start", "stop", "score",
                   "strand", "frame", "attribute", "gene")
# add a column with just gene name to simplify annotation name 
GTF$gene<-gsub(pattern=".*;Name=",replacement= "", x=GTF$gene)
#make a dataframe with bp coverage of all samples
all.bp.coverage <- cbind(as.data.frame(P30$scaffold), P30$position, P30$coverage,
                         P31$coverage, P32$coverage, P33$coverage, P34$coverage,
                         P35$coverage, P36$coverage)
#rename columns
colnames(all.bp.coverage) <- c("scaffold", "position", "P30", "P31",
                               "P32", "P33","P34", "P35" ," P36")
##############################
### Normalizing data 
##############################
#normalize by average coverage across whole genome
# function to determine average coverage of raw bp counts
average.coverage.raw <- function(INPUT,OUTPUT){
  col <- as.numeric(INPUT)
  col.mean<- mean(col, na.rm=TRUE)
  assign(paste(OUTPUT),col.mean, envir = .GlobalEnv)
}
LIST <- as.list(all.bp.coverage[,3:9], all.names=TRUE) # make list with coverage of all samples
  # make sure each element in list is name of sample 
all.coverage.average <- lapply(LIST, 
                FUN  = average.coverage.raw, OUTPUT="cov") # apply function
#calculate mean basepair coverage of all samples
all.coverage.average  <- t(as.data.frame(all.coverage.average))
all.coverage.average <- as.data.frame(all.coverage.average)
mean.all <- mean(all.coverage.average[,1])

# functin to normalize each sample by average coverage
normalize.coverage <- function(FILE, NORM.FACTOR, AVERAGE.FACTOR) {
  temp <- get(x=FILE)
  temp <- data.frame(temp)
  temp <- temp %>%
    mutate(norm.cov = temp[,3] * AVERAGE.FACTOR/ NORM.FACTOR ) 
  assign(FILE,temp, envir = .GlobalEnv)
}
normalize.coverage("P30", all.coverage.average[1,1], mean.all)
normalize.coverage("P31", all.coverage.average[2,1], mean.all)
normalize.coverage("P32", all.coverage.average[3,1], mean.all)
normalize.coverage("P33", all.coverage.average[4,1], mean.all)
normalize.coverage("P34", all.coverage.average[5,1], mean.all)
normalize.coverage("P35", all.coverage.average[6,1], mean.all)
normalize.coverage("P36", all.coverage.average[7,1], mean.all)

write.table(x=P30, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P30.norm.txt",
              quote=FALSE)
write.table(x=P31, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P31.norm.txt",
            quote=FALSE)
write.table(x=P32, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P32.norm.txt",
            quote=FALSE)
write.table(x=P33, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P33.norm.txt",
            quote=FALSE)
write.table(x=P34, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P34.norm.txt",
            quote=FALSE)
write.table(x=P35, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P35.norm.txt",
            quote=FALSE)
write.table(x=P36, file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/P36.norm.txt",
            quote=FALSE)

### normalized by average genome coverage
ALL.COV <- cbind(as.data.frame(P30$scaffold), P30$position, P30$norm.cov,
                 P31$norm.cov, P32$norm.cov, P33$norm.cov, P34$norm.cov,
                 P35$norm.cov, P36$norm.cov)
colnames(ALL.COV) <- c("scaffold", "position", "P30",
                       "P31", "P32", "P33",
                       "P34", "P35", "P36")
ALL.COV$genome.bp <- c(1:nrow(ALL.COV))

# dataframe with RAW non-norm coverage
ALL.COV.RAW <- cbind(as.data.frame(P30$scaffold), P30$position, P30$coverage, P31$coverage, P32$coverage, P33$coverage, P34$coverage, P35$coverage, P36$coverage)
colnames(ALL.COV.RAW) <- c("scaffold", "position", "P30",
                       "P31", "P32", "P33",
                       "P34", "P35", "P36")
#make color panel to color lines for plots below 
SAMPLE.COLORS<- c("P30"="red", "P31" = "red" ,
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
    scale_color_manual(values=SAMPLE.COLORS, name="Sample")+
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


### Calculate the average log2FC difference between pH5, pH7, pH9 groups
calc.log.fc.hk <- function(x) {
  GENE.INFO <- subset(GTF, gene == x)
  GENE.START <- GENE.INFO$start
  GENE.STOP <- GENE.INFO$stop 
  GENE.SCAFFOLD <- as.character(GENE.INFO$scaffold)
  SAMPLE.GENE <- subset(ALL.COV, scaffold ==GENE.SCAFFOLD)
  SAMPLE.GENE <- SAMPLE.GENE %>%
    filter(position > GENE.START & position < GENE.STOP)
  P30_mean <- mean(SAMPLE.GENE$P30)
  P31_mean <- mean(SAMPLE.GENE$P31)
  P32_mean <- mean(SAMPLE.GENE$P32)
  P33_mean <- mean(SAMPLE.GENE$P33)
  P34_mean <- mean(SAMPLE.GENE$P34)
  P35_mean <- mean(SAMPLE.GENE$P35)
  P36_mean <- mean(SAMPLE.GENE$P36)
  pH5_mean <- mean(c(P30_mean, P31_mean))
  pH7_mean <- mean(c(P32_mean, P33_mean, P34_mean))
  pH9_mean <- mean(c(P35_mean, P36_mean))
  pH57_fc <- pH5_mean / pH7_mean
  pH79_fc <- pH9_mean / pH7_mean
  pH59_fc <- pH5_mean / pH9_mean
  temp <- cbind(pH57_fc, pH79_fc, pH59_fc)
  temp <- cbind(temp, pH5_mean, pH7_mean, pH9_mean)
  name1 <- paste("FC",x,sep="_")
  assign(name1,temp,envir = .GlobalEnv)
  }

LIST.HK <- list(GTF[969,10], GTF[2727,10], 
                GTF[3924,10], GTF[3572,10],
                GTF[3864,10])
HK.FC<-sapply(LIST.HK, FUN = calc.log.fc.hk,USE.NAMES=FALSE)
HK.FC.t <- t(HK.FC)
log2HK <- log2(HK.FC.t)
pH57_hk_mean <- mean(HK.FC.t[,1])
pH79_hk_mean <- mean(HK.FC.t[,2])
pH59_hk_mean <- mean(HK.FC.t[,3])

# Prepare data for "operon search"
# calculate average normalized coverage per condition
ALL.COV <- ALL.COV %>%
  mutate(pH5 = P30+P31/2)
ALL.COV <- ALL.COV %>%
  mutate(pH7 = (P32+P33+P34)/3)
ALL.COV <- ALL.COV %>%
  mutate(pH9 = (P35+P36)/2)

ALL.COV <- ALL.COV %>%
  mutate(log2fc_ph57 = log2(pH5/pH7) )
ALL.COV <- ALL.COV %>%
  mutate(log2fc_ph97 = log2(pH9/pH7) )
ALL.COV <- ALL.COV%>%
  mutate(log2fc_ph59 = log2(pH5/pH9) )

pH5_pH7 <- cbind(data.frame(ALL.COV$scaffold, ALL.COV$position))
pH5_pH7 <- cbind(pH5_pH7, data.frame(ALL.COV$log2fc_ph57))
pH5_pH7 <- cbind(pH5_pH7, data.frame(ALL.COV$genome.bp))

pH9_pH7 <- cbind(data.frame(ALL.COV$scaffold, ALL.COV$position))
pH9_pH7 <- cbind(pH9_pH7, data.frame(ALL.COV$log2fc_ph97))
pH9_pH7 <- cbind(pH9_pH7, data.frame(ALL.COV$genome.bp))

# make plot of bp position with logFC <-1 or > 1 
#replace inf with NA
is.na(pH5_pH7) <- do.call(cbind,lapply(pH5_pH7, is.infinite))
#remove NA
pH5_pH7 <- na.omit(pH5_pH7)
pH5_pH7.filtered <- pH5_pH7 %>%
  filter(ALL.COV.log2fc_ph57 > 1 | ALL.COV.log2fc_ph57 < -1)

is.na(pH9_pH7) <- do.call(cbind,lapply(pH9_pH7, is.infinite))
#remove NA
pH9_pH7 <- na.omit(pH9_pH7)
pH9_pH7.filtered <- pH9_pH7 %>%
  filter(ALL.COV.log2fc_ph97 > 1 | ALL.COV.log2fc_ph97 < -1)

color2 <- c(pH5="red", pH9="blue")


plot <- ggplot()
plot + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #geom_point(data=pH5_pH7.unfiltered, aes(x=ALL.COV.bp, y=ALL.COV.log2fc_ph57, color="pH5low"),  size=0.5) +
  #geom_point(data=pH9_pH7.unfiltered, aes(x=ALL.COV.bp, y=ALL.COV.log2fc_ph97, color="pH9low"), size=0.3)+
  ylab("Log2FC") +
  xlab("Genome position")+
  scale_color_manual(values =color2, name = "Condition")+
  ggtitle("Up or down regulated bases compared to pH 7")+
  geom_point(data=pH5_pH7.filtered, aes(x=ALL.COV.bp, y=ALL.COV.log2fc_ph57, color="pH5"),  size=0.5) +
  geom_point(data=pH9_pH7.filtered, aes(x=ALL.COV.bp, y=ALL.COV.log2fc_ph97, color="pH9"), size=0.3) 
  
#########################
#########################
### OPERON SEARCH
########################
########################

### "Sig Operons" Function
##### function takes input data frame in the following format:
  # col 1 = scaffold #
  # col 2 = scaffold position
  # col 3 = log2FC between two conditions 
  # col 4 = bp # 

pH5_pH7.down <- filter(ALL.COV, pH5 >1 & pH7 >1)
pH5_pH7.down <- filter(pH5_pH7.down, log2fc_ph57 < -1)
pH5_pH7.down.f<- cbind(data.frame(pH5_pH7.down$scaffold, pH5_pH7.down$position))
pH5_pH7.down.f <- cbind(pH5_pH7.down.f, data.frame(pH5_pH7.down$log2fc_ph57))
pH5_pH7.down.f <- cbind(pH5_pH7.down.f, data.frame(pH5_pH7.down$bp))

pH5_pH7.up <- filter(ALL.COV, pH5 >1 & pH7 >1)
pH5_pH7.up <- filter(pH5_pH7.up, log2fc_ph57 > 1)
pH5_pH7.up.f<- cbind(data.frame(pH5_pH7.up$scaffold, pH5_pH7.up$position))
pH5_pH7.up.f <- cbind(pH5_pH7.up.f, data.frame(pH5_pH7.up$log2fc_ph57))
pH5_pH7.up.f <- cbind(pH5_pH7.up.f, data.frame(pH5_pH7.up$bp))

pH9_pH7.up <- filter(ALL.COV, pH9 >1 & pH7 >1)
pH9_pH7.up <- filter(pH9_pH7.up, log2fc_ph97 > 1)
pH9_pH7.up.f<- cbind(data.frame(pH9_pH7.up$scaffold, pH9_pH7.up$position))
pH9_pH7.up.f <- cbind(pH9_pH7.up.f, data.frame(pH9_pH7.up$log2fc_ph97))
pH9_pH7.up.f <- cbind(pH9_pH7.up.f, data.frame(pH9_pH7.up$bp))
##### Next chunk of code is for identifying transcriptonal units that meet user-set criteria:
  ### some minimum length (i.e. 1000 bp) where at least 90% of the bp are upreg or downreg 

#set variables
  # user changes LENGTH to some desired threshold 
   LENGTH = 1000
   # user changes INPUT.DF to data frame containing their data
   OPERON.NO = 0
    TEMP.OPERON <-data.frame()
    TENPERCENT = LENGTH * .10 

#input dataframe change colnames
INPUT.DF <- pH9_pH7.up.f
colnames(INPUT.DF) <- c("scaffold", "scaffold.position", "log2fc", "genome.bp")
    for (i in 1:nrow(INPUT.DF)) {
      START = INPUT.DF[i,4]
      THRESHOLD = INPUT.DF[(i+LENGTH),4]
      THRESHOLD.INDEX = i+LENGTH
      UPPER = START + LENGTH + TENPERCENT
      EXACT = START + LENGTH
      OPERON.NO = OPERON.NO + 1 
       # BARE MINIMUM OPERON LENGTH
      if (THRESHOLD >= EXACT & THRESHOLD <= UPPER) {
          SUBSET <- INPUT.DF[i:THRESHOLD.INDEX,]
          SUBSET$operon.no <- OPERON.NO
          TEMP.OPERON <- rbind(TEMP.OPERON, SUBSET)
          TEMP.OPERON <- TEMP.OPERON[!duplicated(TEMP.OPERON[,c('genome.bp')]),]
      }  
    }
TEMP.OPERON[which(!duplicated(TEMP.OPERON[,c('operon.no')])),"sort"] <- c("NONDUP")
TEMP.OPERON[which(is.na(TEMP.OPERON[,"sort"])),"sort"] <- c("DUP")
TEMP.OPERON[which(!duplicated(TEMP.OPERON[,c('operon.no')])),"sort"] <- c("NONDUP")
TEMP.OPERON$start <- ""
COUNTER = 0
i=1
for (i in i:nrow(TEMP.OPERON)) {
  if (TEMP.OPERON[i,6] == "NONDUP" && TEMP.OPERON[i+1,6] == "DUP") {
    COUNTER = COUNTER + 1
    TEMP.OPERON[i,c("start")] <- COUNTER
  }
}
# reassign TEMP.OPERON To new name, export in txt file format 
TEMP.OPERON <- TEMP.OPERON[,-(5:6)]
pH9_pH7_up_1000bp <- TEMP.OPERON
write.table(x=pH9_pH7_up_1000bp, file="~/GoogleDrive/Stenotrophomonas/data/processed/operons/pH9_pH7_up_1000bp.txt",
            quote=FALSE)

