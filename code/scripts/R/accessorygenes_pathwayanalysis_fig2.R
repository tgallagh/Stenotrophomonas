# 11-08-18
# R code for generating plots looking at gene trends across steno phylogenomics


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

'''
for (i in CATEGORY) {
  #make dataframe for 1 pathway at a time, consists of pathway (rownames) and gene components
  genegroup <- as.data.frame(subset(symbol.convert.seed, Category == paste(i)))
  genegroup <- genegroup[!duplicated(genegroup$Gene.x),]
  #name of genegroup
  genegroupname<-as.character(genegroup[1,8])
  #get values for our 2x2 contingency table
  a <-nrow(subset(CF.matches.seed, Category == paste(i))) #a is the number of CF genes in that PATHWAY 
  b <- nrow(genegroup)-a#b is the number of other genes in genome 
  #(#genes from pathway table - a)
  c <-nrow(CF.matches.seed)-a 
    #c is the number of CF genes NOT in that pathway (#sig genes - a)
  d <-nrow(symbol.convert.seed[!duplicated(symbol.convert.seed$Gene.x),]) 
  #d is the number of other genes in genome not in the pathway 
  #make contingency 2x2 table 
  matrix.stats <- matrix(data=c(a,c,b,d), nrow=2, ncol=2)
  ft<-fisher.test(matrix.stats,alternative = "greater")
  temp1 <- as.data.frame(cbind(ft$p.value, genegroupname))
  SEEDS.STATS <- rbind(SEEDS.STATS, temp1)
}
SEEDS.STATS$V1 <- as.numeric(as.character(SEEDS.STATS$V1))
SEEDS.STATS$adj <- p.adjust(SEEDS.STATS$V1, method="fdr")
subset(SEEDS.STATS, V1 < 0.1) '''

'''FLR19SEED <- FLR19SEED[,-(7:9)]
FLR19SEED$strain <- c("FLR19")

K279aSEED <- read.csv("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/K279a/K279a_PROKKA_SEED.csv")
K279aSEED <- K279aSEED[,-1]
K279aSEED$strain <- c("K279a")

average <- FLR19SEED[,5:6]
CAT.SHORT <- CAT.SHORT[!duplicated(CAT.SHORT),]
K279aSEED <- merge(x=K279aSEED, y=CAT.SHORT, by.x="SEEDCategory", by.y="Category",all.x=T)
FLR19SEED <- FLR19SEED[,-(2:3)]
colnames(FLR19SEED) <- c("PROKKA", "gene.name", "Category","CategoryShortened", "strain")
K279aSEED <- K279aSEED[,-2]
colnames(K279aSEED) <- c("Category", "PROKKA", "gene.name", "strain", "CategoryShortened")

SEED <- rbind(K279aSEED, FLR19SEED)
write.csv(SEED, "/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/genome/K279a_FLR19_SEED.csv")'''

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

ggsave("/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/Fig2.eps",
       width=8, height=10,
       device=cairo_ps())
#write.csv(x=PAMatrix, "/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/SupplementalTables/TableS3_SEEDinfo.csv")

 plot <- ggplot()
plot + geom_violin(data=na.omit(Percent.genes.merge), aes(x=CategoryShortened, y=V2, fill=CF.), size=0.9)+
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

dev.off()
embed_fonts(file="/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/Fig2.eps",
            format="eps2write", 
            outfile = "/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/Fig2.eps")
  

#############
############
## dissimilarity of pan-genome presence absence matrix
############
############
rownames(gene_presence_absence) <- gene_presence_absence$Gene
gene.matrix <- gene_presence_absence[,-(2:14)]
gene.matrix.annotated <-gene.matrix[!grepl(gene.matrix$Gene, pattern = "group"),]
gene.matrix.annotated.num <- as.data.frame(sapply(gene.matrix.annotated, FUN=as.numeric))
rownames(gene.matrix.annotated.num) <- rownames(gene.matrix.annotated)
gene.matrix.annotated.num[gene.matrix.annotated.num==1] <- 0
gene.matrix.annotated.num[!gene.matrix.annotated.num==0] <- 1
gene.matrix.annotated.num <- gene.matrix.annotated.num[,-1]

#include non-annotated genes
gene.matrix.all<- gene_presence_absence[,-(2:14)]
rownames(gene.matrix.all) <- rownames(gene.matrix.all)
gene.matrix.all.cut <- as.data.frame(sapply(gene.matrix.all, FUN=as.numeric))
gene.matrix.all.cut[gene.matrix.all.cut==1] <- 0
gene.matrix.all.cut[!gene.matrix.all.cut==0] <- 1
gene.matrix.all.cut <- gene.matrix.all.cut[,-1]
rownames(gene.matrix.all.cut) <- gene.matrix.all$Gene
gene.matrix.all.cut$Gene <- rownames(gene.matrix.all.cut)

### CCA
CAPordination <-capscale(formula=t(gene.matrix.annotated.num)~CF.,data=strain.info,dist="jaccard")

category.scores <- vegan::scores(CAPordination)
category.scores<-data.frame(category.scores$species)
head(category.scores[order(-category.scores$CAP1),])
plot(CAPordination)

### NMDS pan-genome (annotated genes)
MDSplot <- metaMDS(t(gene.matrix.annotated.num), distance="jaccard",
                          k=2)
met.mds.scores <-MDSplot$species
plot(MDSplot)

## NMDS of accessory genes only 
########## ACCESSORY GENES

core.genes <- subset(gene_presence_absence, No..isolates==74)
core.genes <- as.character(core.genes$Gene)
gene.matrix.all.cut <- gene.matrix.all.cut[gene.matrix.all.cut$Gene %ni% core.genes,]
gene.matrix.all.cut <- gene.matrix.all.cut[,-75]
MDSplot <- metaMDS(t(gene.matrix.all.cut), distance="jaccard",
                   k=2)
plot(MDSplot, display=c("sites", "species"), type="text")

data.scores <-as.data.frame(scores(MDSplot))
data.scores$site <- rownames(data.scores) 
data.scores <- merge(x=data.scores, by.x="site", y=strain.info, by.y="Genome.Name")
require(colorRampPalette)

set3 <- colorRampPalette(brewer.pal('Set3',n=12))

data.scores$Patient <- factor(data.scores$Patient, levels = c("AV",
                                                              "BB",
                                                              "CV",
                                                              "FLR19",
                                                              "FMa", 
                                                              "GC",
                                                              "MC",
                                                              "MS",
                                                              "SanG",
                                                              "TG",
                                                              "ZC",
                                                              "environment", "human"))
plot <- ggplot()
plot + geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, color=Patient, shape=Patient),size=3)+
  scale_shape_manual(values=c(15,16,17,3,4,8,15,16,17,3,4,8,15,16,17)) +
  theme(
    plot.background = element_rect(fill="white"),
    panel.background = element_rect(color="black", fill="white"),
    legend.background = element_rect(fill="white", color="black"),
    legend.key = element_rect(fill="white"),
    legend.title = element_blank()
  )


gene.matrix.all.cut.t <- data.frame(t(gene.matrix.all.cut))
gene.matrix.all.cut.t$Strains <- rownames(gene.matrix.all.cut.t)
gene.matrix.all.cut.t <- merge(x=gene.matrix.all.cut.t, by.x="Strains", y=strain.info,  by.y="Genome.Name")
gene.matrix.all.cut.t$CF. <- as.factor(gene.matrix.all.cut.t$CF.)

anosim(dat=gene.matrix.all.cut.t[,2:11979], grouping=gene.matrix.all.cut.t$Patient, distance="jaccard") 

# can we distinguish by habitat using accessory genes
rownames(gene.matrix.all.cut.t) <- gene.matrix.all.cut.t$Strains
accessorydistance <- vegdist(x=gene.matrix.all.cut.t[,2:11979], method="jaccard")


#pcoa
accessorypcoa<-pcoa(D=accessorydistance)
accessorypcoaaxes <- accessorypcoa$vectors
accessorypcoaaxes<- data.frame(accessorypcoaaxes[,1:2])
accessorypcoaaxes$Strains <- rownames(accessorypcoaaxes)
accessorypcoaaxes <- merge(x=accessorypcoaaxes, y= strain.info, by.x="Strains", by.y="Genome.Name")


plot + geom_point(data=accessorypcoaaxes, aes(x=Axis.1, y=Axis.2, color=habitat, shape=Patient))


##### CORE GENOME ALIGNMENT DISSIMILARITY
core.tree<-read.tree("/Volumes/GoogleDrive/My\ Drive/Stenotrophomonas/data/processed/genome/phylogenetic/newtree_alex/accessory_binary_genes.fa.newick")
core.distances <- cophenetic.phylo(core.tree)
coremds <- metaMDS(core.distances)
plot(coremds, display=c("sites"))
core.distances.merge <- data.frame(core.distances)
core.distances.merge$Strains <- rownames(core.distances.merge)
core.distances.merge <- merge(x=core.distances.merge, y=strain.info, by.x="Strains", by.y="Genome.Name")

anosim(dat=core.distances.merge[2:75], grouping=core.distances.merge$Patient, distance="jaccard") 

core.data.scores <-as.data.frame(scores(coremds))
core.data.scores$site <- rownames(core.data.scores) 

core.data.scores <- merge(x=core.data.scores, by.x="site", y=strain.info, by.y="Genome.Name")

plot + geom_point(data=core.data.scores, aes(x=NMDS1, y=NMDS2, color=Patient))


corepcoa<-pcoa(D=core.distances)
corepcoaaxes <- corepcoa$vectors
corepcoaaxes<- data.frame(corepcoaaxes[,1:2])
corepcoaaxes$Strains <- rownames(corepcoaaxes)
corepcoaaxes <- merge(x=corepcoaaxes, y= strain.info, by.x="Strains", by.y="Genome.Name")
corepcoaaxes$Strains<- gsub(corepcoaaxes$Strains, pattern = "Stenotrophomonas_maltophilia_strain_", replacement="")

plot + geom_point(data=corepcoaaxes, aes(x=Axis.1, y=Axis.2, color=Patient, shape=habitat))+
  scale_shape_manual(values=c(15,16,17,3,4,8,15,16,17,3,4,8,15,16,17))

+
  geom_text(data=corepcoaaxes, aes(x=Axis.1, y=Axis.2))

adonis(formula=core.distances.merge[,2:75]~habitat, data = core.distances.merge) # CF vs non-CF


