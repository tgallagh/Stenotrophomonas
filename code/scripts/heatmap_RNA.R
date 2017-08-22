# 7/30/17
# Heatmap of KEGG annotated genes 
require(xlsx)
library("reshape2")
library("dplyr")
library("ggplot2")
library("gplots")
library(scrime)

# heat map of KEGG annotations 
data <- read.xlsx("~/GoogleDrive/Stenotrophomonas/output/tables/heatmap_kegg.xlsx", 
                sheetIndex = 1)
data.matrix <- data[,7:13]
#scale each row so that the mean is zero 
data.matrix.scaled<-rowScales(data.matrix)
data.cut <- cbind(data.frame(data.matrix.scaled), as.character(data$Order), data$Protein.Family,
                  data$Category)
colnames(data.cut) <- c("Acidic \n 1", 
                        "Acidic \n 2",
                        "Neutral \n 1",
                        "Neutral \n 2",
                        "Neutral \n 3",
                        "Basic \n 1",
                        "Basic \n 2",
                        "Order",
                        "Protein.Family",
                        "Category")

data.cut$Protein.Family <- gsub("Bacterial motility", "Bacterial \n motility", data.cut$Protein.Family)
data.cut$Protein.Family <- gsub("Two-component system", "Two-component \n system", data.cut$Protein.Family)
data.cut$Protein.Family <- gsub("Transcription Factors", "Transcription", data.cut$Protein.Family)
data.cut$Protein.Family <- gsub("Translation factor", "Translation", data.cut$Protein.Family)
data.cut$Protein.Family <- gsub("Ribosome", "Translation", data.cut$Protein.Family)
data.cut$Protein.Family <- gsub("Translations", "Translation", data.cut$Protein.Family)

data.melt <- melt(data.cut)
plot <- ggplot(data=data.melt, aes(x=variable, y=as.factor(Order)))
plot + geom_tile(aes(fill=value))+
  facet_grid(Protein.Family~., scales="free", space="free")+
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
        strip.text.y=element_text(angle=0),
        legend.title=element_text(),
        legend.title.align=0.5,
        panel.spacing=unit(.05,"lines"))

#### heatmap by RAST category 
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

plot <- ggplot(data=subset(data.rast.merged.scaled.melt, !Category == "Uncharacterized"), aes(x=variable, y=as.factor(Order)))
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

