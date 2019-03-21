# 9/20/17
# Heatmap by SEED category 
require(xlsx)
library("reshape2")
library("dplyr")
library("ggplot2")
library("gplots")
library(scrime)

#### heatmap by RAST category of all 200 or so differentially expressed genes 
counts.genes<- read.table("/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/edgeR/gene.counts.edgeR.txt",
                     header=T)
colnames(counts.genes) <- c("genes", "Acidic \n 1", 
                        "Acidic \n 2",
                        "Neutral \n 1",
                        "Neutral \n 2",
                        "Neutral \n 3",
                        "Basic \n 1",
                        "Basic \n 2")
data.rast <- read.xlsx("~/GoogleDrive/Stenotrophomonas/output/tables/edgeRallRAST.xlsx",
                       1)
data.rast.merged <- merge(x=data.rast, y=counts.genes, by.x="RAST", by.y="genes")
data.rast.merged.cut <- data.rast.merged[,-(5:10)]

#scale each row so that the mean is zero 
data.rast.merged.matrix <- as.matrix(data.rast.merged.cut[,11:17])
data.rast.merged.matrix.scaled <-rowScales(data.rast.merged.matrix)
data.rast.merged.scaled <- cbind(data.frame(data.rast.merged.matrix.scaled), as.character(data.rast.merged.cut$RAST), 
                                 data.rast.merged.cut$Category)
colnames(data.rast.merged.scaled) <- c("Acidic \n 1", 
                            "Acidic \n 2",
                            "Neutral \n 1",
                            "Neutral \n 2",
                            "Neutral \n 3",
                            "Basic \n 1",
                            "Basic \n 2",
                            "RAST", "Category")

data.rast.merged.scaled <- data.rast.merged.scaled[order(data.rast.merged.scaled$"Acidic \n 1"),]
data.rast.merged.scaled$Order <- c(1:nrow(data.rast.merged.scaled))

data.rast.merged.scaled.melt <- melt(data.rast.merged.scaled, id.vars=c("RAST", "Category", "Order"))

#### all categories
plot <- ggplot(data=subset(data.rast.merged.scaled.melt, !Category == "Uncharacterized"), aes(x=variable, y=as.factor(Order)))
#tiff(file="~/GoogleDrive/Stenotrophomonas/output/figures/paper/seedheatmap.tiff", width=800, height=1000)

plot + geom_tile(aes(fill=value))+
  facet_grid(Category~., scales="free", space="free")+
  scale_fill_gradientn(colours=c("purple", "orange"), name="Relative \n Expression")+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        strip.text.y=element_text(angle=0,size=8),
        legend.title=element_text(),
        legend.title.align=0.5,
        panel.spacing=unit(.0001,"lines"))
#dev.off()


###########################################################
###########################################################
#### only stat sig categories according to fisher exact test 
############################################################
############################################################
SEED <- read.delim("~/GoogleDrive/Stenotrophomonas/output/tables/reformmatedseedannotations.txt")
# list of sig categories from fisher exact test (only got sig categories from pH 5 )
## only suflur metabolism upreg, all others downreg
FT.RESULTS <- c("Sulfur Metabolism", "Virulence, Disease and Defense",
                "Amino Acids and Derivatives", "Motility and Chemotaxis")
SEED.subset <- subset(SEED, Category =="Sulfur Metabolism" |
                        Category == "Virulence, Disease and Defense" |
                        Category == "Amino Acids and Derivatives" |
                        Category == "Motility and Chemotaxis")
counts.genes.sig <- merge(x=SEED.subset, y=counts.genes, by.x="Features", by.y="genes")
#scale each row so that the mean is zero 
data.rast.merged.matrix.scaled <-rowScales(data.rast.merged.matrix)
data.rast.merged.scaled <- cbind(data.frame(data.rast.merged.matrix.scaled), as.character(data.rast.merged.cut$RAST), 
                                 data.rast.merged.cut$Category)
colnames(data.rast.merged.scaled) <- c("Acidic \n 1", 
                                       "Acidic \n 2",
                                       "Neutral \n 1",
                                       "Neutral \n 2",
                                       "Neutral \n 3",
                                       "Basic \n 1",
                                       "Basic \n 2",
                                       "RAST", "Category")



data.rast.merged.scaled.melt$Category <-gsub("Sulfur Metabolism", "Sulfur \n Metabolism", data.rast.merged.scaled.melt$Category) 
data.rast.merged.scaled.melt$Category <-gsub("Virulence, Disease and Defense", "Virulence, Disease \n and Defense", data.rast.merged.scaled.melt$Category) 
data.rast.merged.scaled.melt$Category <-gsub("Amino Acids and Derivatives", "Amino Acids \n and Derivatives", data.rast.merged.scaled.melt$Category) 
data.rast.merged.scaled.melt$Category <-gsub("Motility and Chemotaxis", "Motility and \n Chemotaxis", data.rast.merged.scaled.melt$Category) 


data.rast.merged.scaled.melt$Category <- factor(data.rast.merged.scaled.melt$Category, levels=c("Sulfur \n Metabolism",
                                                                                                "Virulence, Disease \n and Defense",
                                                                                                "Amino Acids \n and Derivatives",
                                                                                                "Motility and \n Chemotaxis"))
plot <- ggplot(data=subset(data.rast.merged.scaled.melt, Category == "Sulfur \n Metabolism" | Category =="Virulence, Disease \n and Defense" | Category == "Amino Acids \n and Derivatives" | Category == "Motility and \n Chemotaxis"), aes(x=variable, y=as.factor(Order)))
#tiff(file="~/GoogleDrive/Stenotrophomonas/output/figures/paper/seedheatmap.tiff", width=800, height=1000)

                                        
plot + geom_tile(aes(fill=value))+
  facet_grid(Category~., scales="free", space="free")+
  scale_fill_gradientn(colours=c("white", "black"), name="Relative \n Expression")+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(color="black", size=14),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right",
        strip.background = element_rect(fill="white", color="black"),
        panel.background=element_rect(fill="white", color="black"),
       # panel.border=element_rect(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        strip.text.y=element_text(angle=0,size=12),
        legend.title=element_text(size=12),
       legend.text=element_text(size=12),
        legend.title.align=0.5,
        panel.spacing=unit(.01,"lines"))
#dev.off()


