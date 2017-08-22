### Annotation and graphs of "operons"
#### TARA
##### R 3.3.1
###### Updated: 7/19/17
require(dplyr)
require(xlsx)
require(ggplot2)

# Read in GTF file of annotated genes from RAST 
GTF <- read.delim(file = "~/GoogleDrive/Stenotrophomonas/data/processed/6666666.230262.gtf",
                  header=FALSE)
GTF$gene <- GTF$V9
#rename columns 
colnames(GTF) <- c("scaffold", "FIG", "type", "start", "stop", "score",
                   "strand", "frame", "attribute", "gene")

up57 <- read.xlsx(file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/pH5_pH7_up_1000bp.xlsm",
                  1)
down57 <- read.xlsx(file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/pH5_pH7_down_1000bp.xlsm",
                  1)
up97 <- read.xlsx(file = "~/GoogleDrive/Stenotrophomonas/data/processed/operons/pH9_pH7_up_1000bp.xlsm",
                  1)

setwd("~/GoogleDrive/Stenotrophomonas/data/processed/operons/")
sample.names <- c("P30", "P31", "P32", "P33", "P34", "P35", "P36")
for (file in sample.names) {
  temp <- read.delim(paste(file,".norm.txt",sep=""),header=FALSE, sep=" ")
  temp <- temp[-1,]
  temp <- temp[,-1]
  colnames(temp) <- c("scaffold", "position", "coverage", "norm.coverage")
  temp$bp <- seq(1:nrow(temp))
  temp$position <- as.numeric(as.character(temp$position))
  assign(file,temp)
}

#function to get annotation info about regions from GTF table
# the input data frame must have the following layout:
  # col1 = scaffold number, col2 = scaffold.position, col3=log2fc, 
  # col4 = genome.bp, col5= operon number
get.annotation.info <- function(df, number) {
  data=subset(df, operon == paste(number))
  scaffold.position<-data$scaffold.position
  scaffold.position.first <- scaffold.position[1]
  scaffold.position.last <- scaffold.position[length(scaffold.position)]
  operon.length <- abs(scaffold.position.last - scaffold.position.first)
  SCAFFOLD <- as.character(data$scaffold[1])
  subsetgtf <- subset(GTF, scaffold == SCAFFOLD)
  subsetgtf<- subsetgtf %>% filter (start > (scaffold.position.first - (operon.length))
                                    & stop < scaffold.position.last + (operon.length))
  # convert scaffold number to genome position to fit onto the graph 
  CONVERT.FACTOR1 <- subset(P30, scaffold == SCAFFOLD & position == scaffold.position[1])
  CONVERT.FACTOR <- as.numeric(CONVERT.FACTOR1$bp - data$scaffold.position[1])
  subsetgtf <- subsetgtf %>% 
    mutate(genome.start = round(CONVERT.FACTOR + start))
  subsetgtf <- subsetgtf %>% 
    mutate(genome.stop = round(CONVERT.FACTOR + stop))
  assign("data", data, envir =  .GlobalEnv)
  assign("subsetgtf", subsetgtf, envir=.GlobalEnv)
}

get.annotation.info(down57, 5)

cols <- c("pH5"="red","pH7"="black","pH9"="blue")
first <- subsetgtf[1,11]
last <- subsetgtf[2,12]

  
plot <- ggplot()
plot + geom_line(data=filter(P30, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH5"))+
  geom_line(data=filter(P31, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH5"))+
  geom_line(data=filter(P32, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH7"))+
  geom_line(data=filter(P33, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH7")) +
geom_line(data=filter(P34, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH7"))+
  geom_line(data=filter(P35, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH9"))+
  geom_line(data=filter(P36, bp >= first & bp <= last), aes(x=bp, y=norm.coverage, color="pH9"))+
  scale_color_manual(values=cols, name="Condition") +
 geom_segment(data=subsetgtf[1,], aes(x=genome.start, xend=genome.stop, y=-10, yend=-10), color="purple",size=4)+
geom_label(data=subsetgtf[1,], aes(label=gene, x=genome.start,y=-100), size=3, fill="purple")+
  geom_segment(data=subsetgtf[2,], aes(x=genome.start, xend=genome.stop, y=-10, yend=-10), color="green",size=4)+
  geom_label(data=subsetgtf[2,], aes(label=gene, x=genome.start,y=-150), size=3, fill="green") 
+
 
   geom_segment(data=subsetgtf[3,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4)+
  geom_label(data=subsetgtf[3,], aes(label=gene, x=genome.start,y=-200), size=3, fill="purple")

+
  
  
  geom_segment(data=subsetgtf[4,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="green",size=4)+
  geom_label(data=subsetgtf[4,], aes(label=gene, x=genome.start,y=-1800), size=3, fill="green") +
 geom_segment(data=subsetgtf[5,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4)+
  geom_label(data=subsetgtf[5,], aes(label=gene, x=genome.start,y=-2200), size=3, fill="purple")+
geom_segment(data=subsetgtf[6,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="purple",size=4)+
  geom_label(data=subsetgtf[6,], aes(label=gene, x=genome.start,y=-2600), size=3, fill="purple")+
  geom_segment(data=subsetgtf[7,], aes(x=genome.start, xend=genome.stop, y=0, yend=0), color="green",size=4)+
  geom_label(data=subsetgtf[7,], aes(label=gene, x=genome.start,y=-3000), size=3, fill="green")
  

  
  
  

