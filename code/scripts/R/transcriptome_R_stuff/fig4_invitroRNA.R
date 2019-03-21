#### smear plots and coverage plots for in vitro rnaseq
### updated 11-14-18

require(xlsx)
require(ggplot2)
require(reshape2)
require(vegan)
require(extrafont)
require(grid)
require(gridExtra)
############################
###########################
# PLOT A, coverage of CPM

# name of individual samples 
SAMPLE.NAMES <- c("P30", "P31", "P32", "P33", "P34", "P35", "P36") 
# conditions of samples in same order as SAMPLE.NAMES vector
SAMPLE.LABELS<- c("pH5", "pH5", "pH7", "pH7", "pH7", "pH9", "pH9")
#Total number of aligned reads, in same order as SAMPLE.NAMES vector
ALIGNED.READS<-c(3249685, 2411591, 3471746, 2187908, 3607402, 5079716, 4023699)
SAMPLE.INFO <- cbind(SAMPLE.NAMES,SAMPLE.LABELS, ALIGNED.READS)
#set working directory with bedtools coverage file 
WORKING.DIRECTORY= "/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/bowtie2/coverage/"
setwd(WORKING.DIRECTORY)
# for loop for reading in all data files with sample name listed in sample.names
# note, suffix of coverage file name must be "_final.coverage.txt" 
# and prefix must be sample name
for (file in SAMPLE.NAMES) {
  temp <- read.delim(paste(file,"_final.coverage.txt",sep=""),header=FALSE)
  colnames(temp) <- c("scaffold", "position", "coverage")
  temp$Type <- c("in vitro")
  assign(file,temp)
}

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

temp<-c()
for (sample in SAMPLE.NAMES) {
  temp <- get(sample)
  temp$bp <- rownames(temp)
  temp$sample <- sample
  temp.list<-slidingwindow(1000, temp[,3])
  df <- data.frame(temp.list[[1]], temp.list[[2]], temp.list[[3]])
  colname<-c("bp","binned_cpm","sd")
  colnames(df)<-colname
  df2 <- merge(x=df, y=temp, by="bp")
  df.name <- paste(sample,"bins", sep="_")
  assign(df.name,df2)
}
#make a dataframe with bp coverage of all samples
P30_bins$Type <- c("in vitro")
all.bp.coverage.invitro <- cbind(as.data.frame(P30_bins$scaffold), P30_bins$bp, P30_bins$Type, P30_bins$coverage,
                                 P31_bins$coverage, P32_bins$coverage, P33_bins$coverage, P34_bins$coverage,
                                 P35_bins$coverage, P36_bins$coverage)
#rename columns
colnames(all.bp.coverage.invitro) <- c("scaffold", "bp", "Type", "P30", "P31",
                                       "P32", "P33","P34", "P35" ," P36")


rm(P30)
rm(P31)
rm(P32)
rm(P33)
rm(P34)
rm(P35)
rm(P36)

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
LIST <- as.list(all.bp.coverage.invitro[,4:10], all.names=TRUE) # make list with coverage of all samples
# make sure each element in list is name of sample 
all.coverage.average <- lapply(LIST, 
                               FUN  = average.coverage.raw, OUTPUT="cov") # apply function
#calculate mean basepair coverage of all samples
all.coverage.average  <- t(as.data.frame(all.coverage.average))
all.coverage.average <- as.data.frame(all.coverage.average)
mean.all <- mean(all.coverage.average[,1])
# function to normalize each sample by average coverage
normalize.coverage <- function(FILE, NORM.FACTOR, AVERAGE.FACTOR) {
  temp <- get(FILE)
  temp <- data.frame(temp)
  temp <- temp %>%
    mutate(norm.cov = temp[,3] * AVERAGE.FACTOR/ NORM.FACTOR ) 
  assign(FILE,temp, envir = .GlobalEnv)
}
normalize.coverage("P30_bins", all.coverage.average[1,1], mean.all)
normalize.coverage("P31_bins", all.coverage.average[2,1], mean.all)
normalize.coverage("P32_bins", all.coverage.average[3,1], mean.all)
normalize.coverage("P33_bins", all.coverage.average[4,1], mean.all)
normalize.coverage("P34_bins", all.coverage.average[5,1], mean.all)
normalize.coverage("P35_bins", all.coverage.average[6,1], mean.all)
normalize.coverage("P36_bins", all.coverage.average[7,1], mean.all)


### normalized by average genome coverage
ALL.COV.NORM.INVITRO <- cbind(as.data.frame(P30_bins$scaffold), P30_bins$bp, P30_bins$Type,P30_bins$norm.cov,
                              P31_bins$norm.cov, P32_bins$norm.cov, P33_bins$norm.cov, P34_bins$norm.cov,
                              P35_bins$norm.cov, P36_bins$norm.cov)
colnames(ALL.COV.NORM.INVITRO) <- c("scaffold", "bp", "Type", "P30",
                                    "P31", "P32", "P33",
                                    "P34", "P35", "P36")
ALL.COV.NORM.INVITRO$genome.bp <- c(1:nrow(ALL.COV.NORM.INVITRO))
# dataframe with RAW non-norm coverage

ALL.COV.RAW <- cbind(as.data.frame(P30_bins$bp), P30_bins$coverage, P31_bins$coverage, P32_bins$coverage, P33_bins$coverage, P34_bins$coverage, P35_bins$coverage, P36_bins$coverage)
colnames(ALL.COV.RAW) <- c( "position", "Acidic \n 1",
                            "Acidic \n 2", "Neutral \n 1", "Neutral \n 2",
                            "Neutral \n 3", "Basic \n 1", "Basic \n 2")


ALL.COV.RAW.m <- melt(ALL.COV.RAW, id.vars = c("position"))


#make color panel to color lines for plots below 


COLORS <-  c(acidic="#e66101", neutral="#cccccc", basic="#5e3c99",
             unbuffered="black")

SAMPLE.COLORS<- c( "Acidic \n 1" = "#e66101" ,
                   "Acidic \n 2" = "#e66101",
                   "Neutral \n 1" = "#4c4c4c",
                   "Neutral \n 2" = "#4c4c4c",
                   "Neutral \n 3" = "#4c4c4c",
                   "Basic \n 1" = "#5e3c99", 
                   "Basic \n 2" = "#5e3c99") 

ALL.COV.RAW.m$Mbp <- ALL.COV.RAW.m$position/1000000

plotA <- ggplot() +
  geom_line(data=ALL.COV.RAW.m, aes(x=as.numeric(Mbp), y=value, color=variable))+
  scale_color_manual(values=SAMPLE.COLORS, name="Condition") +
  facet_grid(variable~.) +
  theme(axis.text=element_text(color="black", size=16),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid =element_line(color="black"),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(), legend.position="none",
        strip.background=element_rect(fill="white", color="black"),
        strip.text=element_text(size=16))+
  labs(x=expression(paste(italic("S. maltophilia"),"FLR19 genome position (Mbp)"))) + 
  labs(y=c("Coverage of RNAseq reads"))

#highcounts <- subset(ALL.COV.RAW.m, value > 5000)
#highcounts <- merge(x=highcounts, y=P30_bins, by.x="position",
                    by.y="bp")



### Bray Curtis and PCOA of just transcriptomes
invitro.rna <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/HTSeq/steno_htseq.xlsx",1)
rownames(invitro.rna) <- invitro.rna$gene
total.counts <- colSums(invitro.rna[,2:8])
INVITRO.CPM <- t((t(invitro.rna[,2:8]) * 1000000) / total.counts)
INVITRO.CPM <-data.frame(INVITRO.CPM)
INVITRO.CPM <- INVITRO.CPM[1:4378,]

distance.matrix.invitro<-vegdist(t(INVITRO.CPM), method="bray")
distance.matrix.coordinates.invitro<-pcoa(distance.matrix.invitro, rn=NULL)
biplot(distance.matrix.coordinates.invitro)

distance.matrix.coordinates.invitro$Group <- factor(distance.matrix.coordinates.invitro$Group,
                                                         levels=c("acidic", "neutral", "basic"))

distance.matrix.axes.invitro<-distance.matrix.coordinates.invitro$vectors
distance.matrix.coordinates.invitro <- data.frame(cbind(distance.matrix.axes.invitro,
                                                        c("acidic", "acidic",
                                                          "neutral", "neutral", "neutral", 
                                                          "basic", "basic")))

distance.matrix.coordinates.invitro$Samples <- rownames(distance.matrix.coordinates.invitro)
distance.matrix.coordinates.invitro.melt <- melt(distance.matrix.coordinates.invitro, id.vars = c("V7", "Samples"))


plotB<-ggplot() + geom_point(data=distance.matrix.coordinates.invitro,aes(x=as.numeric(as.character(Axis.1)), 
                                                               y=as.numeric(as.character(Axis.2)), fill=V7,
                                                               shape=V7), size=5)+
  scale_fill_manual(values=COLORS, name="Group")+
  scale_shape_manual(values=SHAPES, name="Group")+
  scale_x_continuous(limits = c(-0.15, 0.15))+
 # scale_y_continuous(limits=c(-0.15,0.15))+
  xlab("Axis 1 (80.0%)")+
  ylab("Axis 2 (7.0%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),axis.text=element_text(size=16,colour="black"),
        legend.position=c(0.4,0.9),
        legend.text = element_text(size=12),
        legend.title=element_blank(),
        legend.key = element_rect(colour = "white", fill="white"),
        legend.background = element_rect(color="black", fill="white"))+
  ylim(-0.05, .1)


'''invitroaxis <- read.xlsx(file = "~/GoogleDrive/Stenotrophomonas/output/tables/stenotranscriptome_braycurtisdistance.xlsx", 2)
invitroaxis.m <- melt(invitroaxis)
plot + geom_point(data=subset(invitroaxis.m, variable == "Axis.2"),aes(x=NA., 
                                                                       y=as.numeric(as.character(value)), color=Group), size=4)+
  scale_color_manual(values=COLORS, name="Group")+
  scale_shape_manual(values=SHAPES, name="Group")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=16), axis.title.x=element_blank(), axis.text=element_text(size=16,colour="black"),
        legend.text = element_text(size=16),
        legend.title=element_text(size=16),
        legend.key = element_rect(colour = "white", fill="white"))+
  ylab("Bray-Curtis Axis 1 distance to \n averaged neutral transcriptomes")'''


######### SMEAR PLOTS, edgeR output

all57 <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/edgeR/all.5vs7.genes.xlsm", 1)
all57$Threshold <- factor(all57$Threshold, levels = c("Not sig", "FDR", "FDR and FC"))
plotC <- ggplot() + geom_point(data = all57, aes(x=logCPM, y=logFC, color = Threshold))+
  scale_color_manual(values = c("black", "black", "#e6550d"))+
  scale_shape_manual(values= c(17, 20))+    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(),        
legend.position ="none",
axis.text=element_text(size=16, color="black"),
axis.title=element_text(size=16),
axis.line = element_line(colour = "black"))+
  labs(y=expression(paste(log["2"],FC)))+
  labs(x=expression(paste(log["2"],CPM)))

#all79 <- read.xlsx("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/all.7vs9.genes.xlsx", 1)

all79 <- read.delim("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/edgeR/all.7vs9.genes.txt")

all79$Threshold <- factor(all79$Threshold, levels = c("Not sig", "FDR", "FDR and FC"))
plot<-ggplot()

plotD <- ggplot() + geom_point(data = all79, aes(x=logCPM, y=logFC, color = Threshold))+
  scale_color_manual(values = c("black", "black", "purple"))+
  scale_shape_manual(values= c(17, 20))+    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),        
      legend.position="none",
        axis.text=element_text(size=16, color="black"),
        axis.title=element_text(size=16),
        axis.line = element_line(colour = "black"))+
  labs(y=expression(paste(log["2"],FC)))+
  labs(x=expression(paste(log["2"],CPM)))







#### GROB AND EXPORT

tableA<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=(-6))
tableB<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=(-2))
tableC<-textGrob(c("c"), gp=gpar(fontsize=40), vjust=(-2))
tabled<-textGrob(c("d"), gp=gpar(fontsize=40), vjust=(-2))

#tableblank <- textGrob(c(""))

PLOTRNA <-arrangeGrob(
    arrangeGrob(arrangeGrob(tableA, plotA, ncol=2, widths=c(1,12))),
    arrangeGrob(arrangeGrob(tableB, plotB, ncol=2, widths=c(1,12)),
              arrangeGrob(tableC, plotC, ncol=2, widths=c(1,12)),
               arrangeGrob(tabled, plotD, ncol=2, widths=c(1,12)), ncol=3,nrow=1),
   nrow=2, ncol=1, heights=c(2,1)
)
  


ggsave("/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/final/Fig4.eps", 
       device="eps", PLOTRNA, width=11, height=10)
embed_fonts(file="/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/final/Fig4.eps",
           format="eps2write", 
           options="-dDEVICEWIDTHPOINTS=900 -dDEVICEHEIGHTPOINTS=900",
            outfile = "/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/final/Fig4.eps")



