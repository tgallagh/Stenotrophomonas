TABLE <- read.csv("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/gene_presence_absence.csv")
TABLE <- TABLE[,-2]
TABLE <- TABLE[,-(6:12)]
require(reshape2)
require(dplyr)
require(ggplot2)
require(gridExtra)

TABLE.m <- melt(TABLE, id.vars = colnames(TABLE[1:6]))

setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/CFdata/Whiteley/new_aligned/")
# read in sputum E and sputum F
SputumE <- read.delim(file="SputumE.aligned.sam", header=F)
SputumE <- SputumE[-(1:2),]
SputumE <- as.data.frame(SputumE[,3])
SputumE$Sample <- c("SputumE")
SputumE$ID <- SputumE[,1]
SputumE$ID <- gsub(SputumE$ID, pattern ="|gnl.*", replacement="")
SputumE$ID <- gsub(SputumE$ID, pattern ="\\|", replacement="")
SputumE.merged <- merge(x=SputumE, y=TABLE.m, by.x="ID", by.y="value", all.x=TRUE)
SputumE.merged <- SputumE.merged[,-(7:8)]
SputumE.merged <- SputumE.merged[,-2]
SputumE.merged.hits <- SputumE.merged %>%
  group_by(Annotation, No..isolates, Sample, Gene, Avg.group.size.nuc) %>%
  tally()
rm(SputumE)
rm(SputumE.merged)
write.table(SputumE.merged.hits, "SputumE.Hits.txt", quote=F, sep="\t", row.names = F)                                     





SputumF <- read.delim(file="SputumF.aligned.sam", header=F)
SputumF <- SputumF[-(1:2),]
SputumF <- as.data.frame(SputumF[,3])
SputumF$Sample <- c("SputumF")
SputumF$ID <- SputumF[,1]
SputumF$ID <- gsub(SputumF$ID, pattern ="|gnl.*", replacement="")
SputumF$ID <- gsub(SputumF$ID, pattern ="\\|", replacement="")
SputumF.merged <- merge(x=SputumF, y=TABLE.m, by.x="ID", by.y="value", all.x=TRUE)
SputumF.merged <- SputumF.merged[,-(7:8)]
SputumF.merged <- SputumF.merged[,-2]
SputumF.merged.hits <- SputumF.merged %>%
  group_by(Annotation, No..isolates, Sample, Gene, Avg.group.size.nuc) %>%
  tally()
rm(SputumF)
rm(SputumF.merged)
write.table(SputumF.merged.hits, "SputumF.Hits.txt", quote=F, sep="\t", row.names = F)                                     

## get invitro samples
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/invitro_align_final/")

Acidic1 <- read.delim("P30.aligned.sam", header=F)
Acidic1$Sample <- c("Acidic1")
Acidic1 <- as.data.frame(Acidic1[-(1:2),])
Acidic1$ID <- Acidic1[,1]
Acidic1$ID <- gsub(Acidic1$ID, pattern ="|gnl.*", replacement="")
Acidic1$ID <- gsub(Acidic1$ID, pattern ="\\|", replacement="")
Acidic1.merged <- merge(x=Acidic1, y=TABLE.m, by.x="ID", by.y="value", all.x=TRUE)
Acidic1.merged <- Acidic1.merged[,-(7:8)]
Acidic1.merged.hits <- Acidic1.merged %>%
  group_by(Annotation, No..isolates, Sample, Gene, Avg.group.size.nuc) %>%
  tally()
rm(Acidic1)
rm(Acidic1.merged)
write.table(Acidic1.merged.hits, "Acidic1.Hits.txt", quote=F, sep="\t", row.names = F)                                     





Acidic2 <- read.delim("P31.aligned.sam", header=F)
Acidic2$Sample <- c("Acidic2")

Neutral1 <- read.delim("P32.aligned.sam", header=F)
Neutral1$Sample <- c("Neutral1")

Neutral2 <- read.delim("P33.aligned.sam", header=F)
Neutral2$Sample <- c("Neutral2")
Neutral2 <- as.data.frame(Neutral2[-(1:2),])
Neutral2$ID <- Neutral2[,1]
Neutral2$ID <- gsub(Neutral2$ID, pattern ="|gnl.*", replacement="")
Neutral2$ID <- gsub(Neutral2$ID, pattern ="\\|", replacement="")
Neutral2.merged <- merge(x=Neutral2, y=TABLE.m, by.x="ID", by.y="value", all.x=TRUE)
Neutral2.merged <- Neutral2.merged[,-(7:8)]
Neutral2.merged.hits <- Neutral2.merged %>%
                                       group_by(Annotation, No..isolates, Sample, Gene, Avg.group.size.nuc) %>%
                                       tally()
rm(Neutral2)
rm(Neutral2.merged)
#write.table(Neutral2.merged.hits, "Neutral2.Hits.txt", quote=F, sep="\t", row.names = F)                                     


Neutral3 <- read.delim("P34.aligned.sam", header=F)
Neutral3$Sample <- c("Neutral3")

Basic1 <- read.delim("P35.aligned.sam", header=F)
Basic1$Sample <- c("Basic1")
Basic1$ID <- Basic1[,1]
Basic1$ID <- gsub(Basic1$ID, pattern ="|gnl.*", replacement="")
Basic1$ID <- gsub(Basic1$ID, pattern ="\\|", replacement="")
Basic1.merged <- merge(x=Basic1, y=TABLE.m, by.x="ID", by.y="value", all.x=TRUE)
Basic1.merged <- Basic1.merged[,-(7:8)]
Basic1.merged.hits <- Basic1.merged %>%
  group_by(Annotation, No..isolates, Sample, Gene, Avg.group.size.nuc) %>%
  tally()
write.table(Basic1.merged.hits, "Basic1.Hits.txt", quote=F, sep="\t", row.names = F)                                     


Basic2 <- read.delim("P36.aligned.sam", header=F)
Basic2$Sample <- c("Basic2")
Basic2$ID <- Basic2[,1]
Basic2$ID <- gsub(Basic2$ID, pattern ="|gnl.*", replacement="")
Basic2$ID <- gsub(Basic2$ID, pattern ="\\|", replacement="")
Basic2.merged <- merge(x=Basic2, y=TABLE.m, by.x="ID", by.y="value", all.x=TRUE)
Basic2.merged <- Basic2.merged[,-(7:8)]
Basic2.merged.hits <- Basic2.merged %>%
  group_by(Annotation, No..isolates, Sample, Gene, Avg.group.size.nuc) %>%
  tally()
write.table(Basic2.merged.hits, "Basic2.Hits.txt", quote=F, sep="\t", row.names = F)                                     





ALL.DATA <- rbind(SputumE, SputumF)
rm(SputumE, SputumF)

                                       

 SputumE <- subset(ALL.DATA.merged.cut.hits, Sample == "SputumE")
 SputumE$RPKM <- SputumE$n * 1000000 / 162117
 SputumE$RPKM <-  SputumE$RPKM / SputumE$Avg.group.size.nuc
 SputumE$FC <- SputumE$RPKM / 17.59672
 
 SputumF <- subset(ALL.DATA.merged.cut.hits, Sample == "SputumF")
 SputumF$RPKM <- SputumF$n * 1000000 / 1298001
 SputumF$RPKM <-  SputumF$RPKM / SputumF$Avg.group.size.nuc
 SputumF$FC <- SputumF$RPKM / 4.656811
 
 

require(ggplot2)
plot <- ggplot()
plot+geom_histogram(data=ALL.DATA.merged, aes(variable), stat="count")



ALL.DATA.merged.cut.hits <- ALL.DATA.merged.cut.hits[order(ALL.DATA.merged.cut.hits$No..isolates),]
ALL.DATA.merged.cut.hits$Order <- c(1:nrow(ALL.DATA.merged.cut.hits))

# make table for converting the PROKKA ID to rast ID and make seed table
FLR19 <- TABLE[,-(2:5)]
FLR19 <- FLR19[,1:3]
FLR19 <- FLR19[,-2]

# SEED Catgeories
FLR19.SEEDS <- FLR19
SEED<-read.csv("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/RAST/seedannotation.csv")
FLR19.SEEDS <- merge(x=FLR19.SEEDS, y=SEED, by.x="RAST", by.y="Features", all.x=T)


# SEED categories (FLR19 RAST annotation )
ALL.DATA.merged.cut.hits.seedcategories <- merge(x=ALL.DATA.merged.cut.hits,
                                           y=FLR19.SEEDS, 
                                           by.x="Gene", by.y="Gene")
ALL.DATA.merged.cut.hits.seedcategories<- ALL.DATA.merged.cut.hits.seedcategories[order(ALL.DATA.merged.cut.hits.seedcategories$CategoryShortened, ALL.DATA.merged.cut.hits.seedcategories$No..isolates),]

ALL.DATA.merged.cut.hits.seedcategories$Order<-c(1:nrow(ALL.DATA.merged.cut.hits.seedcategories))
ALL.DATA.merged.cut.hits.seedcategories$CategoryShortened <- as.character(ALL.DATA.merged.cut.hits.seedcategories$CategoryShortened)
ALL.DATA.merged.cut.hits.seedcategories[["CategoryShortened"]][is.na(ALL.DATA.merged.cut.hits.seedcategories[["CategoryShortened"]])] <- c("No Category")



ALL.DATA.merged.cut.hits.seedcategories.subset <- 
  subset(ALL.DATA.merged.cut.hits.seedcategories,
         CategoryShortened == "Amino Acids" |
           CategoryShortened == "Carbohydrates" |
           CategoryShortened == "Cell Wall" |
           CategoryShortened == "Cofactors" | 
           CategoryShortened == "DNA Metabolism" |
           CategoryShortened == "Fatty Acids" |
           CategoryShortened == "Nucleosides and Nucleotides"  |
           CategoryShortened == "Protein Metabolism" |
           CategoryShortened == "Respiration" | 
           CategoryShortened =="RNA Metabolism"|
           CategoryShortened =="Stress Response" |
           CategoryShortened =="Virulence" |
           CategoryShortened == "No Category")

ALL.DATA.merged.cut.hits.seedcategories.subset$CategoryShortened <- 
  factor(ALL.DATA.merged.cut.hits.seedcategories.subset$CategoryShortened, 
      levels= c("Amino Acids"
    , "Carbohydrates"
    , "Cell Wall"
    , "Cofactors"  
      , "DNA Metabolism" 
      , "Fatty Acids" 
      , "Nucleosides and Nucleotides"  
      , "Protein Metabolism" 
      , "Respiration"  
      ,"RNA Metabolism"
      ,"Stress Response" 
      ,"Virulence" 
      , "No Category"
  ))
      
plot <- ggplot(data=ALL.DATA.merged.cut.hits.seedcategories)
heatmap<-
  plot + geom_tile(aes(x=Sample, y=as.character(Order), fill=n))+
  facet_grid(CategoryShortened~., space="free", scale="free", switch="x")+
  scale_fill_gradient("log2 \n number of \n reads", low="#ff6600", high="#660066")+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(),
        axis.title.y=element_blank(),
        legend.position="left",
       panel.background=element_rect(fill="gray", color="black"),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
      #  plot.background=element_rect(fill="black"),
      strip.background = element_rect(color="black", fill="white"),
        strip.text.y=element_text(angle=0),
        legend.title=element_text(),
        legend.title.align=0.5,
        panel.spacing=unit(.05,"lines"))


Genecounts <-  
  plot+
  geom_col(aes(y=No..isolates, x=Order))+ylab(label="Number of genomes")+
  coord_flip()+
facet_grid(CategoryShortened~., space="free", scale="free")+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(),
        axis.title.y=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_rect(color="black", fill="white"),
         plot.background=element_rect(fill="white"),
        strip.text.y=element_blank(),
        legend.title=element_text(),
        legend.title.align=0.5,
        panel.spacing=unit(.05,"lines"))  



grid.arrange(heatmap, Genecounts, nrow=1, ncol=2,widths=c(3,1))




FLR19.lookup <- read.delim("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/FLR19.PROKKA.gff", header=F)
FLR19.lookup <- subset(FLR19.lookup, V3 =="CDS")
FLR19.lookup$FASTA <- paste(FLR19.lookup$V1, FLR19.lookup$V4+1, sep=":")
FLR19.lookup$FASTA <- paste(FLR19.lookup$FASTA, FLR19.lookup$V5, sep="-")
FLR19.lookup$ID <- FLR19.lookup$V9
FLR19.lookup$ID <- gsub(FLR19.lookup$ID, pattern="ID=", replacement="")
FLR19.lookup$ID <- gsub(FLR19.lookup$ID, pattern=";Parent=.*", replacement="")
FLR19.lookup<-FLR19.lookup[,10:11]

ALL.DATA <- merge(x=ALL.DATA, y=FLR19.lookup, by.x="V3", by.y="FASTA", all.x=TRUE)

ALL.DATA$ID <- ALL.DATA$V3

                 

