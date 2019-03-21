### Bray Curtis and PCOA of metabolomics and transcriptomics and metatranscriptomic
#### updated 11/18/17
#####  R 3.3.1

require(vegan)
# vegan to get bray curtis distance matrix
# pcoa on distance matrix with ape 
require(ape)
require(reshape2)
require(ggplot2)
require(dplyr)

COLORS <-  c(acidic="#e66101", neutral="#cccccc", basic="#5e3c99",
             blank="black")

SHAPES<-c(acidic=22, neutral=21, basic=24,
         blank=23)

### Bray Curtis and PCOA of metabolomes 
metabolomics <- read.delim("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/metabolomics/mx251577_Smsub.txt")

distance.matrix.metabolomics<-vegdist(t(metabolomics[,-1]), method="bray")
distance.matrix.coordinates.metabolomics<-pcoa(distance.matrix.metabolomics, rn=NULL)



biplot(distance.matrix.coordinates.metabolomics)
distance.matrix.axes.metabolomics<-distance.matrix.coordinates.metabolomics$vectors
Group2 <- c("acidic","acidic","acidic",
            "neutral", "neutral", "neutral",
            "basic", "basic", "basic",
            "blank", "blank", "blank")
distance.matrix.axes.metabolomics <- data.frame(cbind(distance.matrix.axes.metabolomics, Group2))

distance.matrix.axes.metabolomics$Group2 <- factor(distance.matrix.axes.metabolomics$Group2,
                                                         levels=c("acidic", "neutral", "basic", "blank"))

pcoaplot<-ggplot()
pcoaplot <-pcoaplot + geom_point(data=distance.matrix.axes.metabolomics,aes(x=as.numeric(as.character(Axis.1)), 
                                                               y=as.numeric(as.character(Axis.2)), fill=Group2,
                                                               shape=Group2),color="black", size=5)+
  scale_fill_manual(values=COLORS, name="Group")+
  scale_shape_manual(values=SHAPES, name="Group")+
  scale_x_continuous(limits = c(-0.4, 0.4))+
  scale_y_continuous(limits=c(-0.4,0.4))+
  xlab("Axis 1 (71%)")+
  ylab("Axis 2 (15%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),axis.text=element_text(size=16,colour="black"),
        legend.text = element_text(size=16),
        legend.title=element_text(size=16),
        legend.key = element_rect(colour = "white", fill="white"))


'''### Bray Curtis and PCOA of metabolomes 
### metabolomes without blank 
distance.matrix.metabolomics<-vegdist(t(metabolomics[,2:10]), method="bray")
distance.matrix.coordinates.metabolomics<-pcoa(distance.matrix.metabolomics, rn=NULL)
biplot(distance.matrix.coordinates.metabolomics)
distance.matrix.axes.metabolomics<-distance.matrix.coordinates.metabolomics$vectors
Group2 <- c("acidic","acidic","acidic",
            "neutral", "neutral", "neutral",
            "basic", "basic", "basic",
            "blank", "blank", "blank")
distance.matrix.axes.metabolomics <- data.frame(cbind(distance.matrix.axes.metabolomics, Group2))
plot<-ggplot()
plot + geom_point(data=distance.matrix.axes.metabolomics,aes(x=as.numeric(as.character(Axis.1)), 
                                                             y=as.numeric(as.character(Axis.2)), color=Group2,
                                                             shape=Group2), size=4)+
  scale_color_manual(values=COLORS, name="Group")+
  scale_shape_manul(values=SHAPES, name="Group")+
  scale_x_continuous(limits = c(-0.4, 0.4))+
  scale_y_continuous(limits=c(-0.4,0.4))+
  xlab("Axis 1 (54.3%)")+
  ylab("Axis 2 (23.2%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),axis.text=element_text(size=14,colour="black"),
        legend.text = element_text(size=14),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "white", fill="white")
'''

mx<-read.csv("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/processed/metabolomics/mx251577_Sm_flipped.csv")
data.melted<-melt(mx)
mx2 <- mx[,-1]

#permanova
mx2.permanova<- adonis(formula = mx2[,2:708] ~ mx2[,1], 
                       data =mx2, permutations = 999, method ="bray")
mx2.t<-t(mx2)
mx2.t[1,1:13]<-rownames(mx2.t)
# Pairwise comparisons
library("PMCMR")
pairwise.permanova(dist(mx2[ ,2:708], "euclidian"), mx2[,1], nperm = 999,
                   progress = TRUE, p.method = "fdr")



### kruskal wallis test = non parametric anova 
list_of_mets <- as.character(unique(data.melted$variable))
kruskal_wallis_sum <- data.frame()
for(name in list_of_mets) {
  subsetted <- subset(data.melted, variable == name)
  summary<-kruskal.test(subsetted$value ~ subsetted$X.1)
  summary <- unlist(summary)
  p<-summary["p.value"]
  p <- as.numeric(p)
  kruskal_wallis_sum <- append(x=kruskal_wallis_sum,values=p)
}


kruskal_output<-cbind(as.character(list_of_mets), kruskal_wallis_sum)
kruskal_output <- as.data.frame(kruskal_output)
kruskal_output$kruskal_wallis_sum <- as.numeric(as.character(kruskal_output$kruskal_wallis_sum))
kruskal_output<-mutate(.data=kruskal_output,neglogp=-log(kruskal_wallis_sum,base=10)) 
#-log(.05,base=10)

#anything > 1.30103 is sig
filtered1<- filter(kruskal_output, neglogp >= 1.30103) %>%
  mutate(sig = c("yes"))
filtered2<- filter(kruskal_output, neglogp < 1.30103) %>%
  mutate(sig = c("no"))
kruskal_output<-rbind(filtered1, filtered2)
#kruskal_output<-arrange(kruskal_output, V1)
rowname <- rownames(kruskal_output)
kruskal_output <- cbind(rowname, kruskal_output)

scatterplot <- ggplot(kruskal_output, aes(y=neglogp, x=rowname)) 
scatterplot + geom_point(aes(color=sig))+
  ylab("-log(p-value)")+
  xlab("Metabolites")+
  ggtitle("Kruskal Wallis Test for Sig Metabolites")+
  theme_classic()+
  geom_hline(yintercept=1.30103)



filtered_names<-filtered1$V1
filtered_names<-as.character(filtered_names)
posthoc_sum_df<- vector()
for(name in filtered_names) {
  subsetted <- subset(data.melted, variable == name)
  posthoc_sum<-posthoc.kruskal.dunn.test(x=subsetted$value,g=subsetted$X.1,
                                         p.adjust.method = "BH") 
  posthoc_sum <- unlist(posthoc_sum)
  ph7_ph5<-as.numeric(posthoc_sum["p.value1"])
  ph9_ph5<-as.numeric(posthoc_sum["p.value2"])
  TH_ph5<-as.numeric(posthoc_sum["p.value3"])
  ph9_ph7<-as.numeric(posthoc_sum["p.value5"])
  TH_ph7<-as.numeric(posthoc_sum["p.value6"])
  TH_ph9<-as.numeric(posthoc_sum["p.value9"])
  binded<-cbind(ph7_ph5,ph9_ph5)
  binded <- cbind(binded,TH_ph5)
  binded<-cbind(binded, ph9_ph7)
  binded<-cbind(binded,TH_ph7)
  binded<-cbind(binded,TH_ph9)
  posthoc_sum_df <- rbind(binded,posthoc_sum_df)
}
posthoc_sum_df<-cbind(rev(filtered_names), posthoc_sum_df)



#### average to get log FC

data.melted.sum <- 
  data.melted %>%
  group_by(variable, X.1)%>%
  summarise(mean_int=mean(value))

data.melted.sum.cast <- dcast(data.melted.sum, variable ~ X.1)

data.melted.sum.cast <- data.melted.sum.cast %>%
  mutate(logpH9vspH7 = log((`SmpH9` /`SmpH7`), base=2))

data.melted.sum.cast <- data.melted.sum.cast %>%
  mutate(logFCpH5vspH7= log((`SmpH5` /`SmpH7`), base=2))

data.melted.sum.cast <- data.melted.sum.cast %>%
  mutate(logFCpH5vsTH = log((`SmpH5` /`TH`), base=2))

data.melted.sum.cast <- data.melted.sum.cast %>%
  mutate(logFCpH9vsTH = log((`SmpH9` /`TH`), base=2))

data.melted.sum.cast <- data.melted.sum.cast %>%
  mutate(logFCpH5vsTH = log((`SmpH7` /`TH`), base=2))

#write.xlsx(data.melted.sum.cast, 
           file = "/Users/Tara/GoogleDrive/Sm_metabolomics/data/data.melted.sum.cast.logfc.xlsx")

### merge logFC with p-values

merged<- merge(x=data.melted.sum.cast, y=posthoc_sum_df, by.x="variable", by.y="V1")

rownames(merged) <- merged$variable 
attach(merged)


acidicdata<-merged
basicdata<-merged

###### compare pH 5 to pH 7
acidicdata$ph7_ph5 <- as.numeric(as.character(acidicdata$ph7_ph5))

acidicplot <- ggplot(acidicdata, 
                    aes(y=-log10(ph7_ph5), x=logFCpH5vspH7))
acidicplot<-acidicplot+geom_point()+
  ylab("-log10(p)")+
  xlab("log2FC(Acidic/Neutral)")+
  geom_hline(yintercept=1.30103)+
  theme(axis.text=element_text(color="black", size=20),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid =element_line(color="black"), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="white", color="black"),
        strip.text=element_text(size=13))+
  geom_vline(xintercept=c(-1,1), linetype="dotted")+
  annotate("text", x=-1.5, y=1.65, label="putrescine", size=5)+
  annotate("text", x=1.5, y=1.5, label="X129225", size=5)




###### compare pH 9 to pH 7
merged$ph9_ph7 <- as.numeric(as.character(merged$ph9_ph7))
labels<-subset(basicdata, ph9_ph7 < 0.05 & (logpH9vspH7 <-1))

basicplot <- ggplot(basicdata, 
                 aes(y=-log10(ph9_ph7), x=logpH9vspH7))
basicplot<-basicplot+geom_point()+
  ylab("-log10(p)")+
  xlab("log2FC(Basic/Neutral)")+
  geom_hline(yintercept=1.30103)+
  theme(axis.text=element_text(color="black", size=20),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid =element_line(color="black"), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="white", color="black"),
        strip.text=element_text(size=13))+
  geom_vline(xintercept=c(-1,1), linetype="dotted")+
  annotate("text", x=1.540428, y=2, label="gylcolic acid", size=5)+
  annotate("text", x=3.5, y=1.3, label="isothreonic \n acid ", size=5)+
annotate("text", x=2.4, y=1.5, label="ribonic acid ", size=5) +
  annotate("text", x=1, y=1.5,label="salicylic.acid", size=5 )+
  annotate("text", x=2.171141, y=1.7,label="X17128", size=5 )+
  annotate("text", x=2.7, y=1.35,label="X3200", size=5 )+
  annotate("text", x=-1.2, y=1.5,label="methionine", size=5 )+
  xlim(-1.5,4)


  geom_text(data=labels, aes(x=labels$logpH9vspH7, y=-log10(ph9_ph7), 
                             label=labels$variable), check_overlap = TRUE)



  
  require(grid)
  require(gridExtra)

tablea<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=c(-1,0), hjust=c(.5,0))
tableB<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=c(-2.4,0), hjust=c(.5,0))
tableC<-textGrob(c("c"), gp=gpar(fontsize=40), vjust=c(-2,0), hjust=c(.5,0))
  table <- textGrob(c(""))

  grid.arrange(
  arrangeGrob(tablea, pcoaplot,table, ncol=3, widths=c(1,6,2)),
  arrangeGrob(tableB, acidicplot, ncol=2, widths=c(1,12)),
  arrangeGrob(tableC, basicplot, ncol=2, widths=c(1,12)),
              ncol=1, heights=c(1,2,2))
  
##################
  ################
  

putrescine_plot <- ggplot(subset(data.melted, variable == "putrescine"))
putrescine_plot+ geom_boxplot(aes(x=X.1, y=value))+
  xlab("Group")+
  ylab("putrescine intensity")+
  ggtitle("Boxplot of normalized putrescine intensity")

X134401_plot <- ggplot(subset(data.melted, variable == "X134401"))
X134401_plot+ geom_boxplot(aes(x=NA..1, y=value))+
  xlab("Group")+
  ylab("X134401 intensity")+
  ggtitle("Boxplot of normalized X134401 intensity")


X129225_plot <- ggplot(subset(data.melted, variable == "X129225"))
X129225_plot+ geom_boxplot(aes(x=NA..1, y=value))+
  xlab("Group")+
  ylab("X129225 intensity")+
  ggtitle("Boxplot of normalized X129225 intensity")






################
###############
### individual plots


SIG.Metabolites <- subset(data.melted,
                          variable == "glycolic.acid" | 
                            variable ==  "isothreonic.acid" | 
                        variable == "methionine" | variable == "putrescine"  |
                          variable ==  "ribonic.acid" |
                          variable == "salicylic.acid"  | variable== "X129225" |
                          variable == "X17128" | variable== "X3200")


SIG.Metabolites$X.1 <- gsub(SIG.Metabolites$X.1,pattern="Sm", replacement="")
SIG.Metabolites$X.1 <- gsub(SIG.Metabolites$X.1,pattern="pH5", replacement="pH 5")
SIG.Metabolites$X.1 <- gsub(SIG.Metabolites$X.1,pattern="pH7", replacement="pH 7")
SIG.Metabolites$X.1 <- gsub(SIG.Metabolites$X.1,pattern="pH9", replacement="pH 9")

SIG.Metabolites$variable <- gsub(SIG.Metabolites$variable,pattern="\\.", replacement=" ")


ggsave("/Volumes/GoogleDrive/My\ Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/FigS3.eps", 
       device="eps", height=7, width=9)

ggplot()+
  geom_boxplot(data=SIG.Metabolites, aes(x=X.1, y=value))+
  facet_wrap("variable", scales="free")+
  ylab("Normalized intensity")+
  theme(axis.text.x=element_text(color="black", size=14, family="Arial"),
        axis.text.y=element_text(color="black", size=14,  family="Arial"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14, family="Arial"),
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major =element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="white", color="black"),
        strip.text=element_text(size=14, angle=0,  family="Arial"))

dev.off()

embedFonts(
  file=  "/Volumes/GoogleDrive/My\ Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/FigS3.eps",
  format="eps2write",
  outfile="/Volumes/GoogleDrive/My\ Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/FigS3.eps",
  options = "-dEPSCrop")

