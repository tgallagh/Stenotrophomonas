require(reshape2)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(xlsx)
TABLE <- read.csv("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/gene_presence_absence.csv")
TABLE <- TABLE[,-2]
TABLE <- TABLE[,-(6:12)]
TABLE.m <- melt(TABLE, id.vars = colnames(TABLE[1:6]))

# mtrx alignments
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/CFdata/Whiteley/new_aligned/")
FILES <- list.files()
FILES <- FILES[grep(pattern="Hits", FILES)]
temp<-c()
ALL.DATA <-c()
for (SAMPLE in FILES) {
  temp <- read.delim(SAMPLE)
  lib.size<-sum(temp$n)
  temp$RPKM <- temp$n * 1000000 / lib.size
  temp$RPKM <- temp$RPKM / temp$Avg.group.size.nuc
  temp$Normalized <- temp$RPKM / temp$No..isolates
  recA <- subset(temp, Gene =="recA")
  recA <- as.numeric(recA[8])
  temp$FC <- temp$RPKM / recA
  temp <- temp[,-(1:3)]
  temp <- temp[-(2:4)]
  temp <- temp[,-2]
  colnames(temp) <- c("Gene", paste(SAMPLE))
  assign(paste(SAMPLE), temp)
}

# invitro samples
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/invitro_align_final/")
FILES <- list.files()
FILES <- FILES[grep(pattern="Hits", FILES)]

temp<-c()
ALL.DATA <-c()
for (SAMPLE in FILES) {
  temp <- read.delim(SAMPLE)
  lib.size<-sum(temp$n)
  temp$RPKM <- temp$n * 1000000 / lib.size
  temp$RPKM <- temp$RPKM / temp$Avg.group.size.nuc
  temp$Normalized <- temp$RPKM / temp$No..isolates
  recA <- subset(temp, Gene =="recA")
  recA <- as.numeric(recA[8])
  temp$FC <- temp$RPKM / recA
  temp <- temp[,-(1:3)]
  temp <- temp[-(2:4)]
  temp <- temp[,-2]
  colnames(temp) <- c("Gene", paste(SAMPLE))
  assign(paste(SAMPLE), temp)
}

ALL.DATA.table <- merge(x=Acidic1.Hits.txt, y=Acidic2.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=Neutral1.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=Neutral2.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=Neutral3.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=Basic1.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=Basic2.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=SputumE.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=SputumF.Hits.txt, by="Gene", all=T)
TABLE.cut <- TABLE[,1:6]
ALL.DATA.table <-merge(x=ALL.DATA.table, y=TABLE.cut, all.x=T, by="Gene")
ALL.DATA.table<- ALL.DATA.table[order(ALL.DATA.table$No..isolates),]
ALL.DATA.table$Order<-c(1:nrow(ALL.DATA.table))
ALL.DATA.table <- ALL.DATA.table[-6655,]
ALL.DATA.norm <- melt(ALL.DATA.table, id.vars = c("Annotation", "No..isolates", "Gene", "Avg.group.size.nuc", "Order",
                                                  "Avg.sequences.per.isolate", "No..sequences"))


######################################
######################################
### RPKM
######################################
######################################
### FIGURE 5


#setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/invitro_align_final/")
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/CFdata/Whiteley/new_aligned/")

FILES <- list.files()
FILES <- FILES[grep(pattern="Hits", FILES)]

temp<-c()
for (SAMPLE in FILES) {
  temp <- read.delim(SAMPLE)
  lib.size<-sum(temp$n)
  temp$RPKM <- temp$n * 1000000 / lib.size
  temp$RPKM <- temp$RPKM / temp$Avg.group.size.nuc
  #temp$Normalized <- temp$RPKM / temp$No..isolates
  #recA <- subset(temp, Gene =="recA")
  #recA <- as.numeric(recA[8])
  #temp$FC <- temp$RPKM / recA
  temp <- temp[,-(1:3)]
  temp <- temp[-(2:3)]
  #temp <- temp[,-2]
  colnames(temp) <- c("Gene", paste(SAMPLE))
  assign(paste("RPKM_",SAMPLE, sep=""), temp)
}



ALL.DATA.table <- merge(x=RPKM_Acidic1.Hits.txt, y=RPKM_Acidic2.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_Neutral1.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_Neutral2.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_Neutral3.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_Basic1.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_Basic2.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_SputumE.Hits.txt, by="Gene", all=T)
ALL.DATA.table <- merge(x=ALL.DATA.table, y=RPKM_SputumF.Hits.txt, by="Gene", all=T)
TABLE.cut <- TABLE[,1:6]
ALL.DATA.table <-merge(x=ALL.DATA.table, y=TABLE.cut, all.x=T, by="Gene")
ALL.DATA.table<- ALL.DATA.table[order(ALL.DATA.table$No..isolates),]
ALL.DATA.table$Order<-c(1:nrow(ALL.DATA.table))
ALL.DATA.table <- ALL.DATA.table[-6655,]
ALL.DATA.table.RPKM <- ALL.DATA.table
ALL.DATA.RPKM <- melt(ALL.DATA.table, id.vars = c("Annotation", "No..isolates", "Gene", "Avg.group.size.nuc", "Order",
                                             "Avg.sequences.per.isolate", "No..sequences"))

#View(ALL.DATA.norm)
ALL.DATA.norm$Sample <- ALL.DATA.norm$variable
ALL.DATA.norm$Sample <- gsub(ALL.DATA.norm$Sample, pattern=".Hits.txt", replacement="")

ALL.DATA.RPKM$Sample <- gsub(ALL.DATA.RPKM$variable, patter= ".Hits.txt", replacement="")
plot <- ggplot(data=ALL.DATA.RPKM)

heatmap<-
  plot + geom_tile(aes(x=Sample, y=Order, fill=log2(value)))+
 # facet_grid(CategoryShortened~., space="free", scale="free", switch="x")+
  scale_fill_gradient("log2\n(RPKM)", low="#ff6600", high="#660066")+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(color="black", size=12, angle=45,hjust = c(0.75)),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(color="white"),
        axis.title.y=element_blank(),
        legend.position="left",
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="gray", color="black"),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        #  plot.background=element_rect(fill="black"),
        strip.background = element_rect(color="black", fill="white"),
        strip.text.y=element_text(angle=0),
        legend.title=element_text(size=12),
        legend.title.align=0.5,
        panel.spacing=unit(.05,"lines"))+
  geom_hline(yintercept=c(4497))


Genomes <-  
  plot+
  geom_col(aes(y=as.numeric(as.character(No..isolates))/9, x=Order))+
  ylab(label="Number of \n genomes \n")+
  coord_flip()+
 #facet_grid(CategoryShortened~., space="free", scale="free")+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(color="black", size=14),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(size=14),
        axis.title.y=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_rect(color="black", fill="white"),
        plot.background=element_rect(fill="white"),
        strip.text.y=element_blank(),
        legend.title=element_text(),
        legend.title.align=0.5,
        panel.spacing=unit(.05,"lines"))  +
  geom_vline(xintercept=c(4497))



grid.arrange(heatmap, Genomes, nrow=1, ncol=2,widths=c(4,1))


#### PCOA
require(ape)
require(vegan)
### Bray-Curtis and Distance Matrix of Transcriptomes and Metatranscriptomes
#distance matrix on the RPKM values
rownames(ALL.DATA.table.RPKM) <- ALL.DATA.table$Gene
ALL.DATA.table.matrix <- ALL.DATA.table.RPKM[,2:10]

#ALL.DATA.table.matrix <- ALL.DATA.table.matrix[,-8]
# replace NA with zero in table
#ALL.DATA.table.matrix[is.na(ALL.DATA.table.matrix)] <- 0
ALL.DATA.table.matrix <- na.omit(ALL.DATA.table.matrix)
colnames(ALL.DATA.table.matrix) <- c("Acidic1", "Acidic2",
                                     "Neutral1", "Neutral2", "Neutral3", "Basic1",
                                     "Basic2",
                                     "SputumE", "SputumF")
Group <- c("acidic", "acidic", "neutral", "neutral", "neutral",
           "basic", "basic",
           "sputum", "sputum")

distance.matrix.fc<-vegdist(t(ALL.DATA.table.matrix), method="euclidean")
distance.matrix.coordinates<-pcoa(distance.matrix.fc)
biplot(distance.matrix.coordinates)
distance.matrix.coordinates<-distance.matrix.coordinates$vectors
distance.matrix.coordinates <- as.data.frame(cbind(distance.matrix.coordinates, Group))
#distance.matrix.coordinates$Group <- factor(distance.matrix.coordinates$Group,
                                                  #  levels=c("Acidic", "Neutral", "Basic", "Sputum"))

COLORS <-  c(acidic="#e66101", neutral="#cccccc", basic="#5e3c99",
            sputum="green")
SHAPES<-c(acidic=22, neutral=21, basic=24,
          sputum=23)


pcoaplot <- ggplot()
pcoaplot<-pcoaplot + geom_point(data=distance.matrix.coordinates,aes(x=as.numeric(as.character(Axis.1)), 
                                                      y=as.numeric(as.character(Axis.2)), fill=Group,
                                                               shape=Group), color="black", size=5)+
  scale_fill_manual(values=COLORS)+
  scale_shape_manual(values=SHAPES)+
  xlab("Axis 1 (84%)")+
  ylab("Axis 2 (9.8%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),axis.text=element_text(size=16,colour="black"),
        legend.text = element_text(size=16),
        legend.title=element_text(size=16),
        legend.key = element_rect(colour = "white", fill="white"))



SAMPLES <- colnames(ALL.DATA.table.matrix)
Group <- c("invitro", "invitro", "invitro", "invitro", 
           "invitro", "invitro", "invitro", "sputum", "sputum")
sample.info <- data.frame(cbind(SAMPLES, Group))
sample.info$Group <- as.factor(sample.info$Group)
sample.info$Condition <- c("acidic", "acidic",
                           "neutral", "neutral", "neutral",
                           "basic", "basic",
                           "sputum", "sputum")
sample.info$pH <- c("acidic", "acidic",
                    "non-acidic", "non-acidic",
                    "non-acidic", "non-acidic", "non-acidic",
                    "acidic", "acidic"
                    )

#https://r-forge.r-project.org/forum/forum.php?thread_id=2758&forum_id=194&group_id=68
# https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata#238962

adonis(distance.matrix.fc ~ Group/Condition,data=sample.info, permutations=999)
Col1 <- c("Source (pH experiment or sputum)", "Condition (acidic, basic, or \n neutral pH, or sputum)","Residuals")
Col2 <- c(1,2,5)
Col3 <- c(56.7 , 5.1, "")
Col4 <- c(0.79, 0.14, 0.07)
Col5 <- c(.002, .064, "")
permanovatable <- data.frame(cbind(Col1, Col2, Col3, Col4, Col5))
colnames(permanovatable) <- c("Distance Matrix ~ Source/Condition", "DF", "F Model", "Variance Explained (R^2)", "Pr(>F")

permanovatable<-grid.arrange(tableGrob(permanovatable,rows=NULL))

require(randomForest)

tableA<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=c(-4,0), hjust=c(.5,0))
tableB<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=c(-2.4,0), hjust=c(.5,0))
tableC<-textGrob(c(""), gp=gpar(fontsize=40), vjust=c(-2,0), hjust=c(.5,0))
tableD<-textGrob(c(""), gp=gpar(fontsize=40), vjust=c(-2,0), hjust=c(.5,0))

FIG5<-arrangeGrob(
  arrangeGrob(tableA, arrangeGrob(heatmap, Genomes, nrow=1, ncol=2,widths=c(4,1)), ncol=2, widths=c(1,12)),
  arrangeGrob(arrangeGrob(tableC,tableB, pcoaplot, tableC, ncol=4, widths=c(1,1,12,2))), nrow=2, heights=c(2,1))

ggsave("/Volumes/GoogleDrive/My\ Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/Fig5.eps", 
       device="eps", FIG5)
dev.off()

embed_fonts(file="/Volumes/GoogleDrive/My\ Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/Fig5.eps",
            format="eps2write", 
            outfile = "/Volumes/GoogleDrive/My\ Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/Fig5.eps")



### Proportion of reads that align to categories
# mtrx alignments
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/CFdata/Whiteley/new_aligned/")
FILES <- list.files()
FILES <- FILES[grep(pattern="Hits", FILES)]
temp<-c()

for (SAMPLE in FILES) {
  temp <- read.delim(SAMPLE)
 temp <- temp[,-(1:3)]
 temp<-temp[,-2]
  colnames(temp) <- c("Gene", paste(SAMPLE))
  assign(paste("counts",SAMPLE, sep="_"), temp)
}

setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/invitro_align_final/")
FILES <- list.files()
FILES <- FILES[grep(pattern="Hits", FILES)]
temp<-c()
for (SAMPLE in FILES) {
  temp <- read.delim(SAMPLE)
  temp <- temp[,-(1:3)]
  temp<-temp[,-2]
  colnames(temp) <- c("Gene", paste(SAMPLE))
  assign(paste("counts",SAMPLE, sep="_"), temp)
}

# make table for converting the PROKKA ID to rast ID and make seed table
FLR19 <- TABLE[,-(2:5)]
FLR19 <- FLR19[,1:3]
FLR19 <- FLR19[,-2]

FLR19.SEEDS <- FLR19
SEED<-read.csv("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/RAST/seedannotation.csv")
RAST.table <- read.csv("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/prokka_RAST_matchCOMPLETE.csv")
RAST.table <- RAST.table[,2:3]
RAST.table$RAST <- gsub(RAST.table$gene.name, pattern="ID=", replacement="")
RAST.table$RAST <- gsub(RAST.table$RAST, pattern=";Name=.*", replacement="")
FLR19.SEEDS <- merge(x=FLR19.SEEDS, y=RAST.table, by.x="FLR19", by.y="V4")

FLR19.SEEDS <- merge(x=FLR19.SEEDS, y=SEED, by.x="RAST", by.y="Features", all.x=T)
FLR19.SEEDS <- FLR19.SEEDS[,-(1:2)]
FLR19.SEEDS <- FLR19.SEEDS[,-(2:3)]
FLR19.SEEDS <- FLR19.SEEDS[,1:2]

TABLE <- merge(x=ALL.DATA.table.RPKM, y=FLR19.SEEDS, by.x="Gene", by.y="Gene", all.x=T)
edgeR <- read.delim("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/output/tables/edgeRallRAST.txt")
edgeR <- edgeR[,-(1:3)]
edgeR <- edgeR[,-(2:8)]
edgeR <- edgeR[,1:2]
edgeR <- edgeR[unique(edgeR$RAST),]
edgeR$RAST <- as.character(edgeR$RAST)
TABLE <- merge(x=TABLE, y =edgeR,
               all.x=T, by="RAST")


LIST <- ls()
LIST <- LIST[grep(LIST, pattern="counts_")]

TEMP <- c()
ALL.PERCENT <- c()
for (SAMPLEDF in LIST) {
  TEMP <- get(SAMPLEDF)
  TEMP <- merge(x=TEMP, y=FLR19.SEEDS, by="Gene", all.y=T)
  TEMP <- na.omit(TEMP)
  TEMP.total <- sum(TEMP[,2])
  TEMP$percent <- TEMP[,2] / TEMP.total
  TEMP$Sample <- paste(SAMPLEDF)
  colnames(TEMP) <- c("Gene", "Counts", "CategoryShortened", "Percent", "Sample")
 TEMP <-  TEMP %>%
    group_by(CategoryShortened) %>%
   summarise(Total.Percent = sum(Percent)) 
 TEMP <- TEMP[order(TEMP$Total.Percent),]
 TEMP$Order <- 1:nrow(TEMP)
 TEMP$Sample <- paste(SAMPLEDF)
 ALL.PERCENT <- rbind(ALL.PERCENT, TEMP)
}



ALL.PERCENT$Sample <- gsub(ALL.PERCENT$Sample, pattern="counts_", replacement="")
ALL.PERCENT$Sample <- gsub(ALL.PERCENT$Sample, pattern=".Hits.txt", replacement="")

ALL.PERCENT$Group <- ALL.PERCENT$Sample
ALL.PERCENT$Group <- gsub(ALL.PERCENT$Group, pattern="[0-9]", replacement="")
ALL.PERCENT$Group <- gsub(ALL.PERCENT$Group, pattern="E", replacement="")
ALL.PERCENT$Group <- gsub(ALL.PERCENT$Group, pattern="F", replacement="")

ALL.PERCENT.table <- dcast(ALL.PERCENT, CategoryShortened~Sample, value.var = "Total.Percent")
rownames(ALL.PERCENT.table) <- ALL.PERCENT.table[,1]
distance.matrix.percent<-vegdist(t(ALL.PERCENT.table[,-1]), method="euclidean")
#distance.matrix.coordinates<-pcoa(distance.matrix.percent)
#biplot(distance.matrix.coordinates)
#adonis(distance.matrix.percent ~ Group + Condition, data=sample.info, permutations=999)
library(ggdendro)
hc<-hclust(distance.matrix.percent, method="complete")
hcplot <- plot(hc)
ggdendrogram(hc, rotate = F, size = 10)
#dhc <- as.dendrogram(hc)
#ddata <- dendro_data(dhc, type="rectangle")
#adonis(distance.matrix.percent ~ Group + Condition, data=sample.info, permutations=999)


'''pca.coodrinates<-prcomp(t(ALL.PERCENT.table[,-1]))
summary(pca.coodrinates)
loadings(pca.coodrinates)
biplot(pca.coodrinates, col=c("black", "blue"))
loadings(pca.coodrinates)
'''

ALL.percent.mean <- ALL.PERCENT %>%
  group_by(Group, CategoryShortened) %>%
  summarise(Average=mean(Total.Percent))


sample.list<-sample.info$SAMPLES
sample.list <- gsub(sample.list, pattern=".Hits.txt", replacement= "")
library(RColorBrewer)
n <- 26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[1:26]


for (i in sample.list) {
  TEMP.DATA <- subset(ALL.PERCENT, Sample==paste(i))
  tempplot<- ggplot(data=TEMP.DATA)
  tempplot <- tempplot + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
   coord_flip()+
    xlab("Category")+
    theme(
      axis.title.y = element_blank(),   legend.position = "none",
      axis.title.x =element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.x=element_blank()
    )+ ylab("Percent of reads")+
    scale_fill_manual(values=col_vector)
  assign(paste(i,".plot", sep=""), tempplot)
}

dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type="rectangle")
p <- ggplot(segment(ddata)) 
p<-p+ geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))
p<- p+ geom_label(data=label(ddata), aes(x=x,y=y, label=label), label.size=NA,size=5)+
  theme(axis.text = element_blank(),
        axis.title=element_blank( ),
        axis.ticks=element_blank(),
        plot.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white"))

dendro <- as.dendrogram(hc)
ddata <- dendro_data(dendro, type="rectangle")
ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  theme_dendro()

Basic2PLOT <- ggplot(data=subset(ALL.PERCENT, Sample=="Basic2"))
Basic2PLOT<-Basic2PLOT + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.y = element_blank(),   legend.position = "none",
    axis.title.x =element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_blank(),
    plot.background = element_rect(fill="white"),
    panel.background = element_rect(fill="white", color="black")
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("Basic2")



SputumEPLOT <- ggplot(data=subset(ALL.PERCENT, Sample=="SputumE"))
SputumEPLOT<-SputumEPLOT + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.y = element_blank(),   legend.position = "none",
    axis.title.x =element_text(size=10),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12),
    plot.background = element_rect(fill="white"),
    panel.background = element_rect(fill="white", color="black")
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("SputumE")





g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend.plot  <- ggplot(data=ALL.PERCENT)
legend.plot <- legend.plot + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
  xlab("Category")+
  theme(legend.pos="right",
        legend.text=element_text(size=12),
    axis.title.y = element_blank(), 
    axis.title.x =element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank())+
  ylab("Percent of reads")+
  scale_fill_manual(values=col_vector, name="Category")+
  guides(fill=guide_legend(ncol=1))

mylegend<-g_legend(legend.plot)


grid.arrange(p, arrangeGrob(Basic2.plot, Basic1.plot,
                            Neutral3.plot,
                            Neutral1.plot,
                            Neutral2.plot,
                            SputumF.plot,
                            Acidic2.plot,
                            Acidic1.plot,
                            SputumE.plot, nrow=9),mylegend,ncol=3)

grid.arrange(p, arrangeGrob(Basic2PLOT, Basic1PLOT, nrow=9),legend.plot,ncol=3)


grid.arrange(arrangeGrob(SputumE.plot, Acidic1.plot, Acidic2.plot,
                         SputumF.plot, Neutral2.plot, Neutral1.plot, Neutral3.plot,
                         Basic1.plot, Basic2.plot))




#### AVERAGES

SputumEPLOT <- ggplot(data=subset(ALL.PERCENT, Sample=="SputumE"))
SputumEPLOT <- SputumEPLOT + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.y = element_blank(),   legend.position = "none",
    axis.title.x =element_text(size=10),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12)
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("SputumE")

Acidicplot<- ggplot(data=subset(ALL.percent.mean, Group=="Acidic"))
Acidicplot <- 
  Acidicplot + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Average), y=Average))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.y = element_blank(),   legend.position = "none",
    axis.title.x =element_blank(),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12)
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("Acidic")

SputumFPLOT <- ggplot(data=subset(ALL.PERCENT, Sample=="SputumF"))
SputumFPLOT <- SputumFPLOT + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.y = element_blank(),   legend.position = "none",
    axis.title.x =element_text(size=10),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12)
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("SputumF")

Neutralplot<- ggplot(data=subset(ALL.percent.mean, Group=="Neutral"))
Neutralplot <- 
  Neutralplot + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Average), y=Average))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.y = element_blank(),   legend.position = "none",
    axis.title.x =element_blank(),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12)
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("Neutral")


Basicplot<- ggplot(data=subset(ALL.percent.mean, Group=="Basic"))
Basicplot <- 
  Basicplot + geom_col(aes(fill=CategoryShortened,x=reorder(CategoryShortened, Average), y=Average))+
  coord_flip()+
  xlab("Category")+
  theme(
    axis.title.x = element_blank(),   legend.position = "none",
    axis.title.y =element_blank(),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12)
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=col_vector)+
  ggtitle("Basic")


dendrogram <-ggdendrogram(hc, rotate = T, size = 8, leaf_labels = TRUE)+
  theme(axis.text=element_text(color="black",size=12))

require(grid)
tablea<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=c(-6,0), hjust=c(.5,0))
tableb<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=c(-6,0), hjust=c(.5,0))


grid.arrange(
  arrangeGrob(tablea,arrangeGrob(Acidicplot,Basicplot,Neutralplot,
              SputumEPLOT,SputumFPLOT, 
                          ncol=3, nrow=2), ncol=2, widths=c(1,12)),
      arrangeGrob(tableb,dendrogram, ncol=2, widths=c(1,12)),
             nrow=2, heights=c(4,1))
             
COLORS <-  c(Acidic="#e66101", Neutral="#cccccc", Basic="#5e3c99",
             Sputum="#a1d99b")


relativeplot<- ggplot(data=subset(ALL.PERCENT))
relativeplot<- relativeplot + geom_bar(stat="identity", position="fill",
                aes(fill=Group,x=reorder(CategoryShortened, Total.Percent), y=Total.Percent))+
  coord_flip()+
  xlab("Category")+
  theme(
    legend.position="left",
    legend.text=element_text(size=12),
    axis.title.x = element_blank(),   
    axis.title.y =element_blank(),
    axis.text.y = element_text(color="black", size=12),
    axis.ticks.x = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=12),
    plot.title=element_text(hjust=0.5, size=12)
  )+ ylab("Percent of reads")+
  scale_fill_manual(values=COLORS)

grid.arrange(arrangeGrob(tablea, relativeplot, ncol=2, widths=c(1,12)),
             arrangeGrob(tableb, dendrogram, ncol=2, widths=c(1,12)),
             ncol=2, widths=c(3,1))

### ACCESSORY GENES
CF.accessory <- read.xlsx("/Volumes/GoogleDrive/My Drive//Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/CF.Genes.xlsx",1)
CF.accessory <- CF.accessory[,1:5]

ALL.DATA.accessory <- merge(x=ALL.DATA.table, y=CF.accessory, by="Gene")
