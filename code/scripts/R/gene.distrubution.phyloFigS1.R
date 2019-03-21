### updated 12-12-18

require(xlsx)
require(ggplot2)
require(reshape2)
require(dplyr)

CF.gene.matrix <- read.xlsx("~/GoogleDrive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/CF.Genes.xlsx",1)
all.gene.matrix<-read.csv("~/GoogleDrive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/gene_presence_absence.csv")
strain.info <- read.csv("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/strain_metadata.csv")
#strain.info <- strain.info[,c(2,32,33)]
strain.info$Genome.Name <- as.character(strain.info$Genome.Name)
strain.info$Patient <- as.character(strain.info$Patient)

strain.info <- rbind(strain.info, c("FLR19", "CF", "CF", "FLR19"))
strain.info$CF. <- as.numeric(strain.info$CF.)
strain.info[is.na(strain.info)] <- 0
strain.info$Genome.Name <- gsub(strain.info$Genome.Name, pattern =" ", replacement="_")
strain.info <- strain.info[strain.info$Genome.Name %in% colnames(gene.matrix.annotated),]
strain.info <- rbind(strain.info, c("Stenotrophomonas_maltophilia_strain_CBF101","environment", "1", "environment"))
strain.info <- rbind(strain.info, c("Stenotrophomonas_sp_KAs_53_strain_KAs_53","environment", "1", "environment"))
######################
## look at antibiotics
######################

all.gene.matrix.antibiotics <- all.gene.matrix[,-(5:14)]

all.gene.matrix.antibiotics<-all.gene.matrix.antibiotics[grepl(all.gene.matrix.antibiotics$Annotation,
                                  ignore.case=T, pattern="antibiotic|resistance|drug"),]

all.gene.matrix.antibiotics.melt <- melt(all.gene.matrix.antibiotics,
                                         id.vars = colnames(all.gene.matrix.antibiotics[,1:5]))
#all.gene.matrix.antibiotics.melt$value[all.gene.matrix.antibiotics.melt$value==""]<-"0"
#all.gene.matrix.antibiotics.melt$value[all.gene.matrix.antibiotics.melt$value!=""]<-"1"

all.gene.matrix.antibiotics.melt <- all.gene.matrix.antibiotics.melt[-which(all.gene.matrix.antibiotics.melt$value == ""), ]


all.gene.matrix.antibiotics.melt$value<-as.numeric(as.character((all.gene.matrix.antibiotics.melt$value)))


all.gene.m<- all.gene.matrix[,-(2:14)]
rownames(all.gene.m) <- all.gene.m[,1]
all.gene.m <- all.gene.m[,-1]
all.gene.m[all.gene.m==""]  <- NA
all.gene.m<-replace(data.frame(lapply(all.gene.m, as.character), stringsAsFactors = FALSE),
          !is.na(all.gene.m), "1")
all.gene.m[is.na(all.gene.m)] <- 0
all.gene.melt <- melt(all.gene.m,measure.vars  = colnames(all.gene.m))

## get average number of genes across 
all.gene.melt.omit <- subset(all.gene.melt, !value == "0")
genes.per.strain <- all.gene.melt.omit %>% group_by(variable) %>% tally()
genes.per.strain$variable <- gsub(genes.per.strain$variable, pattern ="Stenotrophomonas_maltophilia_strain_", replacement="")
genes.per.strain$variable <- gsub(genes.per.strain$variable, pattern ="Stenotrophomonas_maltophilia_", replacement="")
genes.per.strain$variable <- gsub(genes.per.strain$variable, pattern ="Stenotrophomonas_sp_KAs_53_strain_", replacement="")

# get meta-data about the strains
metadata <- read.delim("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/strain_metadata.txt")
metadata$Genome.Name <- gsub(metadata$Genome.Name, pattern ="Stenotrophomonas maltophilia strain ", replacement="")
metadata$Genome.Name <- gsub(metadata$Genome.Name, pattern ="Stenotrophomonas maltophilia ", replacement="")
metadata$Genome.Name <- gsub(metadata$Genome.Name, pattern ="Stenotrophomonas sp. KAs 5-3 strain ", replacement="")
metadata <- rbind(metadata, c("FLR19", "CF"))
## add CBF101 and KAs_53 manually
metadata<- rbind(metadata, c("CBF101", "environment"))
metadata <- rbind(metadata, c("KAs_53", "environment"))

genes.per.strain.merge <- merge(x=genes.per.strain, y=metadata, by.x="variable", by.y="Genome.Name")

avg.genes.per.strain <- mean(genes.per.strain$n)

genes.per.strain.merge$variable <- gsub(genes.per.strain.merge$variable, pattern="_", replacement=" ")
genes.per.strain.merge$habitat <- gsub(genes.per.strain.merge$habitat, pattern="human", replacement="non-CF \n human")

plot1 <- ggplot() + geom_col(data=genes.per.strain.merge, 
      aes(x=reorder(variable, desc(variable)), y=n))+
 geom_hline(yintercept=4285)+
  facet_grid(habitat~., space="free", scale="free")+
  coord_flip()+
  ylab("Number of genes per genome")+
  theme(axis.text.x=element_text(color="black", size=10, family="Arial"),
        axis.text.y=element_text(color="black", size=8,  family="Arial"),
        axis.title.x=element_text(size=14,  family="Arial"),
        axis.title.y=element_blank(),
      panel.grid.minor=element_blank(),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major =element_blank(),
      axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="white", color="black"),
        strip.text.y=element_text(size=12, angle=0,  family="Arial"))
  #ylim(3900,5000)

#############
###### FIG S1B
#############

#### Enrichment of genes in pathways

require(xlsx)
require(dplyr)
require(vegan)
require(ape)
require(ggplot2)
require(extrafont)
require(reshape2)

setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/CFaccessorygenome/")
accessory.genes <- read.xlsx("annotatedgenes.xlsx", 1)
accessory.genes$Gene.1 <- as.character(accessory.genes$Gene.1 )
accessory.genes$Gene.1 <- gsub(accessory.genes$Gene.1, pattern=" ", replacement="")
CFcopies<-accessory.genes[duplicated(accessory.genes$Gene.1),]
CFcopies<-accessory.genes[accessory.genes$Gene.1 %in% CFcopies$Gene.1,]
'%ni%' <- Negate('%in%')
unique <- accessory.genes[accessory.genes$Gene.1 %ni% CFcopies$Gene.1,]

all.genes <- read.csv("allgenes.csv") # this is all genes with some annotation in pan-genome
all.genes.list <- as.character(unique(all.genes$Gene.1)) 

#############
#### SEED ANALYSIS
FLR19SEED<- read.csv("StenoSEEDCategories.csv")
#FLR19SEED <- na.omit(FLR19SEED)

CATEGORY <- as.character(unique(SEED$Category))
SEEDS.STATS <- c()
gene_presence_absence <- read.csv("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/gene_presence_absence.csv")

SEED <- read.csv("../K279a_FLR19_SEEDed.csv")

k279a_flr19_matrix <- gene_presence_absence[,c(1,15,16)]
symbol.convert.seed <- merge(x=all.genes, y=k279a_flr19_matrix, by="Gene")
symbol.convert.seed <- melt(symbol.convert.seed, measure.vars = c("FLR19", "Stenotrophomonas_maltophilia_K279a"))

symbol.convert.seed <- merge(x=symbol.convert.seed, y=SEED,by.x="value", by.y="PROKKA")
symbol.convert.seed$both <- paste(symbol.convert.seed$Gene, symbol.convert.seed$CategoryShortened, sep="_")
symbol.convert.seed <- symbol.convert.seed[!duplicated(symbol.convert.seed$both),]

CF.matches.seed <- merge(x=accessory.genes, y=FLR19SEED , by.x="Gene.1", by.y="Symbol")

PAMatrix <- gene_presence_absence[,-(2:14)]
PAMatrix$Symbol <- gsub(PAMatrix$Gene, pattern="_.*", replacement="")
PAMatrix <- merge(x=symbol.convert.seed, y=PAMatrix, by.x="Gene.1", by.y="Symbol")
#PAMatrix <- PAMatrix[!duplicated(PAMatrix)]
Strains <- colnames(PAMatrix)
Strains <- Strains[12:85]
CAT.SHORT <- FLR19SEED[,5:6]
temp<-c()
Percent.genes <- c()
# get percent of gnes that each strain has within a pathway
for (i in CATEGORY){
  CAT.DS <- paste(i)
  for (i in Strains) {
    STRAIN.NAME <- paste(i)
    Cat.subset <- subset(PAMatrix, Category == paste(CAT.DS))
    Cat.subset <- Cat.subset[!duplicated(Cat.subset$Gene.x),]
    total.counts <- nrow(Cat.subset)
    STRAIN.SUBSET <- Cat.subset[STRAIN.NAME]
    STRAIN.SUBSET[STRAIN.SUBSET==""] <-NA
    STRAIN.SUBSET <- na.omit(STRAIN.SUBSET)
    strain.counts <- nrow(STRAIN.SUBSET)
    prop <- strain.counts/total.counts
    temp <- c(STRAIN.NAME, prop, CAT.DS)
    Percent.genes <- rbind(Percent.genes, temp)
  }
}


strain.info <- read.csv("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/strain_metadata.csv")
#strain.info <- strain.info[,c(2,32,33)]
strain.info$Genome.Name <- as.character(strain.info$Genome.Name)
strain.info$Patient <- as.character(strain.info$Patient)


strain.info <- rbind(strain.info, c("FLR19", "CF", "CF", "FLR19"))
strain.info$CF. <- as.numeric(strain.info$CF.)
strain.info[is.na(strain.info)] <- 0
strain.info$Genome.Name <- gsub(strain.info$Genome.Name, pattern =" ", replacement="_")
strain.info <- strain.info[strain.info$Genome.Name %in% colnames(gene.matrix.annotated),]
strain.info <- rbind(strain.info, c("Stenotrophomonas_maltophilia_strain_CBF101","environment", "1", "environment"))
strain.info <- rbind(strain.info, c("Stenotrophomonas_sp_KAs_53_strain_KAs_53","environment", "1", "environment"))

strain.info$Patient <- gsub(strain.info$Patient, pattern = "Fma", replacement="FMa")

Percent.genes <- as.data.frame(Percent.genes)
Percent.genes$V2 <- as.numeric(as.character(Percent.genes$V2))

Percent.genes$V1 <- gsub(Percent.genes$V1, pattern = "FLR19.y", replacement="FLR19")

Percent.genes.merge <- merge(x=Percent.genes, y=strain.info,
                             by.x="V1", by.y="Genome.Name", all.x=T)

Percent.genes.merge$CF. <- as.character(Percent.genes.merge$CF.)
Percent.genes.merge$CF. <- gsub(Percent.genes.merge$CF., pattern="1", replacement="non-CF")
Percent.genes.merge$CF. <- gsub(Percent.genes.merge$CF., pattern="2", replacement="CF")
Percent.genes.merge$CF. <-factor(Percent.genes.merge$CF., levels = c("non-CF", "CF"))
Percent.genes.merge <- merge(x=Percent.genes.merge, y=CAT.SHORT, by.y="Category", by.x="V3")

temp <-c()
CATEGORY.TTEST <- c()
CATEGORY.SUB <- CATEGORY[-10]
for (i in CATEGORY.SUB) {
  CAT <- i
  pathway1 <- subset(Percent.genes.merge, V3 ==CAT)
  CF.pathway1 <- subset(pathway1, CF.=="CF")
  nonCF.pathway1 <- subset(pathway1, CF.=="non-CF")
  CAT.t<-t.test(x=CF.pathway1$V2, y=nonCF.pathway1$V2)
  print(CAT.t)
  CAT.p <- CAT.t$p.value
  CAT.mean<-data.frame(CAT.t$estimate)
  temp <- cbind(CAT, CAT.p, t(CAT.mean))
  CATEGORY.TTEST <- rbind(CATEGORY.TTEST, temp)
}

CATEGORY.TTEST <- data.frame(CATEGORY.TTEST)
CATEGORY.TTEST$CAT.p <- as.numeric(as.character(CATEGORY.TTEST$CAT.p))
CATEGORY.TTEST$padj <- p.adjust(CATEGORY.TTEST$CAT.p, method="bonf") 

# two categoires enriched for in CF = protein metabolism, and virulence
category.pvalues <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/category_pvalues.xlsx",1)

enrichment<-ggplot() + geom_violin(data=na.omit(Percent.genes.merge), aes(x=CategoryShortened, y=V2, fill=CF.), size=0.9)+
  #geom_point(data=Percent.genes.merge, aes(x=V3, y=V2, color=CF.))+
  facet_grid(CategoryShortened~., scales="free", space="free")+
  coord_flip()+
  theme(axis.line=element_blank(),
        axis.text.x=element_text(color="black", size=12, hjust = c(0.75)),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        legend.position="left",
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        #  plot.background=element_rect(fill="black"),
        strip.background = element_rect(color="black", fill="white"),
        strip.text.y=element_text(angle=0, size=10),
        # panel.spacing=unit(.05,"lines")
  )+
  ylab("Proportion of genes in category") +
  geom_text(data=category.pvalues, show.legend=F,aes(x=CategoryShortened,y=V2), label="*", size=8)+
  scale_fill_manual(values=c("#cccccc", "green"))


################# 
## make panel
################
require(grid)
require(gridExtra)
require(extrafont)
require(extrafontdb)

tableA<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=c(-10,0), hjust=c(.5,0))
tableB<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=c(-10,0), hjust=c(.5,0))

setwd("~/Desktop")

PLOT<-arrangeGrob(
  arrangeGrob(tableA, plot1, ncol=2, widths=c(1,12)),
  arrangeGrob(tableB, enrichment, ncol=2, widths=c(1,12)),
  nrow=1, ncol=2, widths=c(2,3))

ggsave("./FigS1.eps", 
       device="eps", PLOT, width=12, height=9.3)

embedFonts(
  file=  "./FigS1.eps",
  format="eps2write",
  outfile="./FigS1.eps",
  options = "-dEPSCrop")





