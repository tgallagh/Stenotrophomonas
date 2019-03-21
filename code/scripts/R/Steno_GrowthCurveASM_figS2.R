##### updated 11-14-18
#####  R 3.3.1

require(xlsx)
require(reshape2)
require(ggplot2)
require(dplyr)
require(extrafont)
require(grid) 
require(gridExtra)


setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/")

COLORS <-  c(acidic="#e66101", neutral="#cccccc", basic="#5e3c99",
             unbuffered="black")


SHAPES<-c(acidic=22, neutral=24, basic=21,
          unbuffered=23)



#### FIGURE S2E = ASM PLOT

asm <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/asmdata.xlsx",1)
asm.sig <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/asm_significance.xlsx",1)
asm <- na.omit(asm)
asm <- subset(asm, !Experiment==1)
asm.growth.averages <- asm %>%
  group_by(Media, Timepoint) %>%
  summarise(mean=mean(as.numeric(as.character(Concentration))))
asm.growth.averages <-na.omit(asm.growth.averages)

### fig. S2e
asmplot <- ggplot()
asmplot <- asmplot + geom_point(data=subset(asm, !Timepoint==0 & !Experiment==1),
              aes(x=as.numeric(as.character(Timepoint)), y=as.numeric(as.character(Concentration)),
              fill=Media,
  shape=Media),size=2)+
  scale_fill_manual(values=COLORS, name="")+
  scale_shape_manual(values=SHAPES, name="")+
  #facet_grid(Strain~.,scales="free_y")+
  ylab("Concentration \n (CFU/ml)")+
  xlab("Hours")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14,colour="black"),
        panel.background=element_rect(color="black", fill="white"),
        strip.background = element_rect(color="black", fill="white"),
        strip.text.y=element_text(size=12),
        legend.position = "none",
        legend.text=element_text(size=12, color="black"),
        legend.background = element_rect(fill="white", color="white"),
        legend.key = element_rect(fill="white", color="white"))+
  geom_line(data=subset(asm.growth.averages, !Timepoint==0), size=1, aes(x=as.numeric(as.character(Timepoint)),y=mean,
                                                            group=Media, color=Media))+
 scale_color_manual(values=COLORS, name="")+
  scale_y_continuous(trans="log10") +
geom_text(data=asm.sig, show.legend=F,aes(x=x,y=y, color=Media), label="*", size=6)+
  ggtitle("FLR19 growth in ASM")

data <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/pHonly_allstrains.xlsx",1)
data$Condition <- factor(data$Condition,levels=c("acidic", "neutral", "basic",
                                                 "unbuffered", "citric acid", "lactic acid", "sulfuric acid"))
data.growth <- data[,-(7:9)]
data.growth$Strain <- as.factor(gsub(data.growth$Strain, pattern = "Fma 2012", replacement= "FMa 2012"))

data.growth$Strain <- factor(data.growth$Strain, levels =
                               c(
                                 "FLR19", 
                                 "CV 2008", 
                                 "FMa 2012" , 
                                 "GC 2011", 
                                 "ZC 2005",
                                 "ZC 2006",
                                 "K279a",
                                 "NCTC 10257"
                               ))
data.growth.averages <- data.growth %>%
  group_by(Condition, Timepoint, Strain, GRID) %>%
            summarise(mean=mean(as.numeric(as.character(OD))))


#### ACIDS PLOT
  
  ### Fig S2B

  acidsplot <- ggplot()
   acidsplot <- acidsplot+ geom_point(data=data.growth, aes(x=Timepoint, y=OD, fill=Condition,shape=Condition),size=3)+
    #facet_wrap("Strain", scales="free_x")+
    scale_fill_manual(values=c("#e66101", "#cccccc", "#5e3c99","black", "#2E8B57", "#556B2F",
                                "#00FA9A"), name="")+
    facet_grid(Strain~., scales = "free")+
    scale_shape_manual(values=c(22,24,21,18,25,25,25), name="")+
    ylab(expression(OD[600]))+
    xlab("Incubation (h)")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=16),
          axis.text=element_text(size=12,colour="black"),
          panel.background=element_rect(color="black", fill="white"),
          strip.background = element_rect(color="black", fill="white"),
          strip.text.y=element_text(size=10),
          legend.text=element_text(size=12, color="black"),
          legend.background = element_rect(fill="white", color="white"),
          legend.key = element_rect(fill="white", color="white"))+
    scale_color_manual(values=c("#e66101", "#cccccc", "#5e3c99","black", "#2E8B57", "#556B2F",
                                "#00FA9A"), name="")+
    geom_line(data=data.growth.averages, size=0.75, aes(x=as.numeric(as.character(Timepoint)),y=mean,
                                                               group=Condition, color=Condition))
    
  
  ######## pH cross over
  ##### Fig S2C
  
  crossover<-read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/pHcrossover.xlsx",1)
  crossover$Strain <- factor(crossover$Strain, levels =
                                 c(
                                   "FLR19", 
                                   "CV 2008", 
                                   "FMa 2012" , 
                                   "GC 2011", 
                                   "ZC 2005",
                                   "ZC 2006",
                                   "K279a",
                                   "NCTC 10257"
                                 ))
  crossover.averages <- crossover %>%
    group_by(Condition, Timepoint, Strain, Beginning, Crossover) %>%
    summarise(mean=mean(as.numeric(as.character(OD))))
  
  crossplot <- ggplot()
  crossplot <- crossplot + geom_point(data=crossover,
    aes(x=as.numeric(as.character(Timepoint)), y=as.numeric(as.character(OD)), fill=Beginning, shape=Beginning),color="black",size=2)+
    scale_fill_manual(values=COLORS, name="Beginning pH")+
    scale_shape_manual(values=SHAPES, name="Beginning pH")+
    facet_grid(Strain~Crossover,scales="free_y")+
    ylab(expression(OD[600]))+
    xlab("Hours")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=16),
          axis.text=element_text(size=12,colour="black"),
          panel.background=element_rect(color="black", fill="white"),
          strip.background = element_rect(color="black", fill="white"),
          strip.text.y=element_text(size=10),
          legend.text=element_text(size=12, color="black"),
          legend.background = element_rect(fill="white", color="white"),
          legend.key = element_rect(fill="white", color="white"))+
    geom_line(data=na.omit(crossover.averages), size=0.75, aes(x=as.numeric(as.character(Timepoint)),y=mean,
                                                              group=Beginning, color=Beginning))+
    scale_color_manual(values=COLORS, name="Beginning pH")
    #geom_text(data=sig.ph.only, show.legend=F,aes(x=x,y=y, color=Condition), label="*", size=6)
  
  
### pH PLOT
  
data.pH <- data[,-(8:10)]
data.pH <- subset(data.pH, Strain=="FLR19")
data.pH <- subset(data.pH, !Timepoint==48)
data.pH <- subset(data.pH, Condition == "acidic" |
                    Condition == "neutral" |
                    Condition == "basic" |
                    Condition == "unbuffered"
                  )

data.pH <- data.pH[!is.na(data.pH$pH),]

data.pH.averages <- data.pH %>%
  group_by(Condition, Timepoint) %>%
  summarise(mean=mean(as.numeric(as.character(pH))))


pHplot<-ggplot()
pHplot <- pHplot + geom_point(data=data.pH,aes(x=as.numeric(as.character(Timepoint)), y=as.numeric(as.character(pH)), fill=Condition,
                                       shape=Condition),color="black",size=3)+
  scale_fill_manual(values=c("#e66101", "#cccccc", "#5e3c99","black", "#2E8B57", "#556B2F",
                             "#00FA9A"), name="")+
  scale_shape_manual(values=c(22,24,21,18,25,25,25), name="")+
  ylab("pH")+
  xlab("Hours")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14,colour="black"),
        panel.background=element_rect(color="black", fill="white"),
        strip.background = element_rect(color="black", fill="white"),
        strip.text.y=element_text(size=12),
        legend.text=element_text(size=12, color="black"),
        legend.position="none",
        legend.background = element_rect(fill="white", color="white"),
        legend.key = element_rect(fill="white", color="white"))+
  scale_color_manual(values=COLORS,name="")+
  ggtitle("FLR19 change in pH")+
  geom_line(data=na.omit(data.pH.averages), size=0.75, aes(x=as.numeric(as.character(Timepoint)),y=mean,
                                                             group=Condition, color=Condition))





tableA<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=c(-7)) 
tableB<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=c(-7)) 
tableC<-textGrob(c("c"), gp=gpar(fontsize=40)) 
tableD<-textGrob(c("d"), gp=gpar(fontsize=40)) 
table <- textGrob(c(""), gp=gpar(fontsize=40))


setwd("/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoJBactRebuttal/figures/")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#mylegend <- g_legend(ODplot)

top<-arrangeGrob(arrangeGrob(tableA, acidsplot, ncol=2, widths=c(1,12)),
            arrangeGrob(tableB, crossplot, ncol=2, widths=c(1,12)), ncol=2, widths = c(2,3))

bottom<-arrangeGrob(arrangeGrob(tableC, pHplot,table, ncol=3, widths=c(1,6,4)),
                 arrangeGrob(tableD, asmplot, table, ncol=3, widths=c(1,6,2)), ncol=2, widths=c(1,1))
PLOTOD<-arrangeGrob(top, bottom, ncol=1, nrow=2, heights=c(7,2))

ggsave("FigS2.eps", device="eps", PLOTOD, width=15, height=9.3)
embed_fonts(file="./FigS2.eps",format="eps2write", outfile = "./FigS2embed.eps")

