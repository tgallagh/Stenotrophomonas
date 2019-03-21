##### updated 12-12-18
#####  R 3.3.1

require(xlsx)
require(reshape2)
require(ggplot2)
require(dplyr)
require(extrafont)


setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/")
COLORS <-  c(acidic="#e66101", neutral="#808080", basic="#5e3c99",
             unbuffered="black")

SHAPES<-c(acidic=22, neutral=24, basic=21,
          unbuffered=23)

data <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/pHonly_allstrains.xlsx",1)

data$Condition <- factor(data$Condition,levels=c("acidic", "neutral", "basic",
                                                 "unbuffered"))
data.growth <- data[,-(7:9)]
data.growth$Strain <- as.factor(gsub(data.growth$Strain, pattern = "Fma 2012", replacement= "FMa 2012"))

data.growth.averages <- data.growth %>%
  group_by(Condition, Timepoint, Strain, GRID) %>%
            summarise(mean=mean(as.numeric(as.character(OD))))

#### OD PLOT TO LOOK AT GROWTH CURVES IN PH BUFFERED MEDIA
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/")
sig.ph.only <- read.xlsx("pHonly_significance.xlsx",1)



data.growth$Strain <- factor(data.growth$Strain,
                             levels = c("FLR19",
                                        "CV 2008",
                                        "FMa 2012",
                                          "GC 2011",
                                        "ZC 2005",
                                        "ZC 2006",
                                        "K279a",
                                        "NCTC 10257"))

ODplot <- ggplot()
ODplot <- ODplot + geom_point(data=subset(data.growth, Condition == "acidic" |
                                  Condition == "neutral" | Condition == "basic" | Condition == "unbuffered"),aes(x=as.numeric(as.character(Timepoint)), y=as.numeric(as.character(OD)), fill=Condition,
                                                               shape=Condition),color="black",size=2)+
  scale_fill_manual(values=COLORS, name="")+
  scale_shape_manual(values=SHAPES, name="")+
  facet_grid(Strain~.,scales="free_y")+
  ylab(expression(OD[600]))+
  xlab("Hours")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),
    axis.text=element_text(size=14,colour="black"),
        panel.background=element_rect(color="black", fill="white"),
    strip.background = element_rect(color="black", fill="white"),
    strip.text.y=element_text(size=12),
       legend.text=element_text(size=12, color="black"),
    legend.background = element_rect(fill="white", color="white"),
    legend.key = element_rect(fill="white", color="white"),
    legend.position=c(0.5,1))+
geom_line(data=na.omit(data.growth.averages), size=1, aes(x=as.numeric(as.character(Timepoint)),y=mean,
                                              group=Condition, color=Condition))+
  scale_color_manual(values=COLORS, name="")+
  geom_text(data=sig.ph.only, show.legend=F,aes(x=x,y=y, color=Condition), label="*", size=6)+
  ylim(0,2.1)




#### Stats of OD over time
### COMPARE FLR19 pH 5 to pH 7 

STRAIN.LIST <- levels(data.growth$Strain)
all.p <- c()
all.p.final <-c()

data.growth$log2OD <- log2(data.growth$OD) 
'''
### t - test
for (i in STRAIN.LIST) {
  STRAIN.name <- i
  fourhouracidic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "acidic"| Condition =="neutral") & Timepoint=="4")
  fourhouracidict<-t.test(fourhouracidic$OD ~ fourhouracidic$Condition)
  fourhouracidic.p<-data.frame(fourhouracidict$p.value)
  
  eighthouracidic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "acidic"| Condition =="neutral") & Timepoint=="8")
  eighthouracidict<-t.test(eighthouracidic$OD ~ eighthouracidic$Condition)
  eighthouracidic.p<-data.frame(eighthouracidict$p.value)
  
  twelvehouracidic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "acidic"| Condition =="neutral") & Timepoint=="12")
  twelvehouracidict<-t.test(twelvehouracidic$OD ~ twelvehouracidic$Condition)
  twelvehouracidic.p<-data.frame(twelvehouracidict$p.value)
  
  twentyfourhouracidic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "acidic"| Condition =="neutral") & Timepoint=="24")
  twentyfourhouracidict<-t.test(twentyfourhouracidic$OD ~ twentyfourhouracidic$Condition)
  twentyfourhouracidic.p<-data.frame(twentyfourhouracidict$p.value)
  
  fortyeighthouracidic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "acidic"| Condition =="neutral") & Timepoint=="48")
  fortyeighthouracidict<-t.test(fortyeighthouracidic$OD ~ fortyeighthouracidic$Condition)
  fortyeighthouracidic.p<-data.frame(fortyeighthouracidict$p.value)
  
  fourhourbasic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "basic"| Condition =="neutral") & Timepoint=="4")
  fourhourbasict<-t.test(fourhourbasic$OD ~ fourhourbasic$Condition)
  fourhourbasic.p<-data.frame(fourhourbasict$p.value)
  
  eighthourbasic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "basic"| Condition =="neutral") & Timepoint=="8")
  eighthourbasict<-t.test(eighthourbasic$OD ~ eighthourbasic$Condition)
  eighthourbasic.p<-data.frame(eighthourbasict$p.value)
  
  twelvehourbasic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "basic"| Condition =="neutral") & Timepoint=="12")
  twelvehourbasict<-t.test(twelvehourbasic$OD ~ twelvehourbasic$Condition)
  twelvehourbasic.p<-data.frame(twelvehourbasict$p.value)
  
  twentyfourhourbasic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "basic"| Condition =="neutral") & Timepoint=="24")
  twentyfourhourbasict<-t.test(twentyfourhourbasic$OD ~ twentyfourhourbasic$Condition)
  twentyfourhourbasic.p<-data.frame(twentyfourhourbasict$p.value)
  
  fortyeighthourbasic <- subset(data.growth, Strain == paste(STRAIN.name) & (Condition == "basic"| Condition =="neutral") & Timepoint=="48")
  fortyeighthourbasict<-t.test(fortyeighthourbasic$OD ~ fortyeighthourbasic$Condition)
  fortyeighthourbasic.p<-data.frame(fortyeighthourbasict$p.value)
  
  all.p <- cbind(fourhouracidic.p, eighthouracidic.p, twelvehouracidic.p,
    twentyfourhouracidic.p, fortyeighthouracidic.p,
    fourhourbasic.p, eighthourbasic.p, twelvehourbasic.p,
    twentyfourhourbasic.p, fortyeighthourbasic.p)
  all.p <- data.frame(t(all.p))
  all.p$Comparison <- rownames(all.p)
  all.p$Strain <- STRAIN.name
  
  all.p.final <- rbind(all.p, all.p.final)
  
}'''

'''all.pairwise<-c()
all.pairwise.final<-c()
### anova with post-hoc pairwise comparisons
for (i in STRAIN.LIST) {
  STRAIN.name <- i
  temp <- subset(data.growth, !Condition == "unbuffered")
  
  fourhour <- subset(temp, Strain == paste(STRAIN.name) &Timepoint=="4")
  fourhourpairwise<-pairwise.t.test(fourhour$OD, g= fourhour$Condition, p.adj="none")
  fourhourp<-data.frame(fourhourpairwise$p.value)
  fourhourp$variable <- rownames(fourhourp)
  fourhour.m <- melt(fourhourp)
  fourhour.m <- na.omit(fourhour.m)
  fourhour.m$Timepoint <- c(4)
  
  eighthour <- subset(temp, Strain == paste(STRAIN.name) &Timepoint=="8", p.adj="none")
  eighthourpairwise<-pairwise.t.test(eighthour$OD, g= eighthour$Condition)
  eighthourp<-data.frame(eighthourpairwise$p.value)
  eighthourp$variable <- rownames(eighthourp)
  eighthour.m <- melt(eighthourp)
  eighthour.m <- na.omit(eighthour.m)
  eighthour.m$Timepoint <- c(8)
  
  twelvehour <- subset(temp, Strain == paste(STRAIN.name) &Timepoint=="12", p.adj="none")
  twelvehourpairwise<-pairwise.t.test(twelvehour$OD, g= twelvehour$Condition)
  twelvehourp<-data.frame(twelvehourpairwise$p.value)
  twelvehourp$variable <- rownames(twelvehourp)
  twelvehour.m <- melt(twelvehourp)
  twelvehour.m <- na.omit(twelvehour.m)
  twelvehour.m$Timepoint <- c(12)
  
  twentyfourhour <- subset(temp, Strain == paste(STRAIN.name) &Timepoint=="24", p.adj="none")
  twentyfourhourpairwise<-pairwise.t.test(twentyfourhour$OD, g= twentyfourhour$Condition)
  twentyfourhourp<-data.frame(twentyfourhourpairwise$p.value)
  twentyfourhourp$variable <- rownames(twentyfourhourp)
  twentyfourhour.m <- melt(twentyfourhourp)
  twentyfourhour.m <- na.omit(twentyfourhour.m)
  twentyfourhour.m$Timepoint <- c(24)
  
  fortyeighthour <- subset(temp, Strain == paste(STRAIN.name) &Timepoint=="48", p.adj="none")
  fortyeighthourpairwise<-pairwise.t.test(fortyeighthour$OD, g= fortyeighthour$Condition)
  fortyeighthourp<-data.frame(fortyeighthourpairwise$p.value)
  fortyeighthourp$variable <- rownames(fortyeighthourp)
  fortyeighthour.m <- melt(fortyeighthourp)
  fortyeighthour.m <- na.omit(fortyeighthour.m)
  fortyeighthour.m$Timepoint <- c(48)
  
 all.pairwise <- rbind(fourhour.m,
                       eighthour.m,
                       twelvehour.m,
                       twentyfourhour.m,
                       fortyeighthour.m)

  all.pairwise$Strain <- STRAIN.name
  
  all.pairwise.final <- rbind(all.pairwise, all.pairwise.final)
  
}

all.pairwise.final$adj <- p.adjust(all.pairwise.final$value, "bonf")
all.pairwise.final.filter<- subset(all.pairwise.final,adj < 0.05)
all.pairwise.final.filter$both <- paste(all.pairwise.final.filter$variable,all.pairwise.final.filter$variable.1, sep="_" )
all.pairwise.final.filter <- subset(all.pairwise.final.filter, !both=="basic_acidic")
  
  
SFLR19.5to7 <- subset(data.growth, Strain=="FLR19" & (Condition == "acidic" | Condition == "neutral"))
toutput<-t.test(FLR19.5to7$OD ~ FLR19.5to7$Condition)
'''


#### ANTIBIOTICS PLOTS
setwd("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/")
FLR19 <- read.csv("FLR19_Antibiotics.csv")
italy<-read.csv("italy_strains_antibiotics.csv")
NCTC <- read.csv("NCTC_Antibiotics.csv")
K279a<-read.csv("K279a_Antibiotics.csv")

colnames(italy) <- colnames(FLR19)

antibiotics <- rbind(FLR19,
                         italy, NCTC, K279a)

antibiotics$pH <- gsub(antibiotics$pH, pattern="Acidic", replacement="acidic")
antibiotics$pH <- gsub(antibiotics$pH, pattern="Basic", replacement="basic")
antibiotics$pH <- gsub(antibiotics$pH, pattern="Neutral", replacement="neutral")


antibiotics.sig <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/antibiotics.sig.xlsx",1)


antibiotics$Strain <- gsub(antibiotics$Strain, pattern = "CV2008", replacement="CV 2008")
antibiotics$Strain <- gsub(antibiotics$Strain, pattern = "FMA2012", replacement="FMa 2012")
antibiotics$Strain <- gsub(antibiotics$Strain, pattern = "ZC2005", replacement="ZC 2005")
antibiotics$Strain <- gsub(antibiotics$Strain, pattern = "ZC2006", replacement="ZC 2006")
antibiotics$Strain <- gsub(antibiotics$Strain, pattern = "GC2011", replacement="GC 2011")

antibiotics$Strain <- factor(antibiotics$Strain, 
                             levels = c("FLR19",
                                        "CV 2008",
                                        "FMa 2012",
                                        "GC 2011", 
                                        "ZC 2005",
                                        "ZC 2006",
                                        "K279a",
                                        "NCTC 10257"))


antibiotics<-antibiotics[!is.na(antibiotics$OD),]

data.antbx.averages <- antibiotics%>%
  dplyr::group_by(Antibiotic, Concentration, pH, Strain, Timepoint) %>%
  dplyr::summarise(mean=mean(as.numeric(as.character(OD))))

#########################
#### ANTIBIOTICS OD PLOTS
########################

antibiotics$Concentration <- as.numeric(as.character(antibiotics$Concentration))
data.antbx.averages$Concentration <- as.numeric(as.character(data.antbx.averages$Concentration))
antibiotics.sig$x <- as.numeric(as.character(antibiotics.sig$x))

antibiotics.sig$Strain <- factor(antibiotics.sig$Strain, 
                                 levels = c("FLR19",
                                            "CV 2008",
                                            "FMa 2012",
                                            "GC 2011", 
                                            "ZC 2005",
                                            "ZC 2006",
                                            "K279a",
                                            "NCTC 10257"))


antbxplot <- ggplot()+
  facet_grid(Strain~Antibiotic, scales="free_x")+
  geom_text(data=antibiotics.sig, show.legend=F, aes(x=x,y=y, 
                                            color=pH), label="*", size=7)+
  geom_point(data=subset(antibiotics, Timepoint==24), 
             aes(x=Concentration, y=OD, fill=pH,shape=pH),size=2)+
  scale_fill_manual(values=COLORS, name="")+
  scale_shape_manual(values=SHAPES, name="")+
  xlab("Concentration (mg/L)")+
  ylab(expression(OD[600]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),
        axis.text.x=element_text(size=11,color="black", angle=0),
        axis.text.y=element_text(size=12,colour="black"),
        legend.text = element_text(size=13),
        legend.position="none",
        panel.background=element_rect(color="black", fill="white"),
        strip.background = element_rect(color="black", fill="white"),
        strip.text=element_text( size=12))+
  geom_line(data=subset(data.antbx.averages, Timepoint==24), 
            show.legend=F,size=1,
            aes(x=Concentration,y=mean, 
                group=interaction(pH, Antibiotic), color=pH
            ))+
  scale_color_manual(values=COLORS, name="")+
  scale_y_continuous(expand=c(0,.4))


###########################
######### FOLD-CHANGE PLOTS
###########################

Antibiotics.FC <- subset(antibiotics, Timepoint==24)
Antibiotics.FC<- Antibiotics.FC %>% dplyr::group_by(Antibiotic,pH, Strain, DateExptStarted) %>%
  dplyr::mutate(FC = OD / OD[Concentration == 0])

data.antbx.averages.fc <- Antibiotics.FC%>%
  dplyr::group_by(Antibiotic, Concentration, pH, Strain) %>%
  dplyr::summarise(meanFC=mean(as.numeric(as.character(FC))))

  
#Antibiotics.FC$Concentration <- as.factor(as.character(Antibiotics.FC$Concentration))
#Antibiotics.FC$Concentration <- factor(Antibiotics.FC$Concentration,
                                       levels = c(0, 16, 32, 64,
                                                  128, 200, 256, 400, 800))
#data.antbx.averages.fc$Concentration <- as.factor(as.character(data.antbx.averages.fc$Concentration))
#data.antbx.averages.fc$Concentration <- factor(data.antbx.averages.fc$Concentration,
                                       levels = c(0, 16, 32, 64,
                                                  128, 200, 256, 400, 800))
data.antbx.averages.fc$Concentration <- as.numeric(as.character(data.antbx.averages.fc$Concentration))
Antibiotics.FC$Concentration <- as.numeric(as.character(Antibiotics.FC$Concentration))


fc.sig <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/antibiotics_significance.xlsx",1)

fc.sig$Strain <- factor(fc.sig$Strain,
                        levels =c(
                          "FLR19",
                          "CV 2008",
                          "FMa 2012",
                          "GC 2011", 
                          "ZC 2005",
                          "ZC 2006",
                          "K279a",
                          "NCTC 10257"
                        ))

fc.sig$X <- as.numeric(as.character(fc.sig$X))

Antibiotics.FC$Strain <- factor( Antibiotics.FC$Strain,
                                 levels =c(
                                   "FLR19",
                                   "CV 2008",
                                   "FMa 2012",
                                   "GC 2011", 
                                   "ZC 2005",
                                   "ZC 2006",
                                   "K279a",
                                   "NCTC 10257"
                                 ))

FCplot <- ggplot()+
geom_text(data=fc.sig, show.legend=F, aes(x=x,y=c(1.3), 
                    color=Condition), label="*", size=7)+
  geom_point(data=subset(Antibiotics.FC), 
             aes(x=Concentration, y=FC, fill=pH,shape=pH),size=2)+
  scale_fill_manual(values=COLORS, name="")+
  scale_shape_manual(values=SHAPES, name="")+
  ylab("Fold-change relative to control")+
  xlab("Concentration (mg/L)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=16),
        axis.text.x=element_text(size=11,color="black", angle=0),
        axis.text.y=element_text(size=12,colour="black"),
        legend.text = element_text(size=13),
        legend.position="none",
        panel.background=element_rect(color="black", fill="white"),
        strip.background = element_rect(color="black", fill="white"),
        strip.text=element_text( size=12))+
     geom_line(data=subset(data.antbx.averages.fc), 
            show.legend=F,size=1,
            aes(x=Concentration,y=meanFC, 
                group=interaction(pH, Antibiotic), color=pH
  ))+
  scale_color_manual(values=COLORS, name="")+
 scale_y_continuous(expand=c(0,.5))  +
  facet_grid(Strain~Antibiotic, scales="free")
 
   
'''all.pairwise<-c()
all.pairwise.final<-c()
### anova with post-hoc pairwise comparisons for antibiotics data
STRAIN.LIST <- levels(as.factor(antibiotics$Strain))

for (i in STRAIN.LIST) {
STRAIN.name <- i
temp <- subset(antibiotics, !pH == "unbuffered")
temp <- subset(temp, Timepoint=="24")

Gent0 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "0")
Gent0pairwise<-pairwise.t.test(Gent0$OD, g=Gent0$pH, p.adj="none")
Gent0pairwise<-data.frame(Gent0pairwise$p.value)
Gent0pairwise$variable <- rownames(Gent0pairwise)
Gent0pairwise.m <- melt(Gent0pairwise)
Gent0pairwise.m <- na.omit(Gent0pairwise.m)
Gent0pairwise.m$Antibiotic <- c("Gentamicin")
Gent0pairwise.m$Concentration <- c("0")

Gent200 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "200")
Gent200pairwise<-pairwise.t.test(Gent200$OD, g=Gent200$pH, p.adj="none")
Gent200pairwise<-data.frame(Gent200pairwise$p.value)
Gent200pairwise$variable <- rownames(Gent200pairwise)
Gent200pairwise.m <- melt(Gent200pairwise)
Gent200pairwise.m <- na.omit(Gent200pairwise.m)
Gent200pairwise.m$Antibiotic <- c("Gentamicin")
Gent200pairwise.m$Concentration <- c("200")

Gent400 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "400")
Gent400pairwise<-pairwise.t.test(Gent400$OD, g=Gent400$pH, p.adj="none")
Gent400pairwise<-data.frame(Gent400pairwise$p.value)
Gent400pairwise$variable <- rownames(Gent400pairwise)
Gent400pairwise.m <- melt(Gent400pairwise)
Gent400pairwise.m <- na.omit(Gent400pairwise.m)
Gent400pairwise.m$Antibiotic <- c("Gentamicin")
Gent400pairwise.m$Concentration <- c("400")

Gent800 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "800")
Gent800pairwise<-pairwise.t.test(Gent800$OD, g=Gent800$pH, p.adj="none")
Gent800pairwise<-data.frame(Gent800pairwise$p.value)
Gent800pairwise$variable <- rownames(Gent800pairwise)
Gent800pairwise.m <- melt(Gent800pairwise)
Gent800pairwise.m <- na.omit(Gent800pairwise.m)
Gent800pairwise.m$Antibiotic <- c("Gentamicin")
Gent800pairwise.m$Concentration <- c("800")


Tobra0 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "0")
Tobra0pairwise<-pairwise.t.test(Tobra0$OD, g=Tobra0$pH, p.adj="none")
Tobra0pairwise<-data.frame(Tobra0pairwise$p.value)
Tobra0pairwise$variable <- rownames(Tobra0pairwise)
Tobra0pairwise.m <- melt(Tobra0pairwise)
Tobra0pairwise.m <- na.omit(Tobra0pairwise.m)
Tobra0pairwise.m$Antibiotic <- c("Tobramycin")
Tobra0pairwise.m$Concentration <- c("0")

Tobra16 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "16")
Tobra16pairwise<-pairwise.t.test(Tobra16$OD, g=Tobra16$pH, p.adj="none")
Tobra16pairwise<-data.frame(Tobra16pairwise$p.value)
Tobra16pairwise$variable <- rownames(Tobra16pairwise)
Tobra16pairwise.m <- melt(Tobra16pairwise)
Tobra16pairwise.m <- na.omit(Tobra16pairwise.m)
Tobra16pairwise.m$Antibiotic <- c("Tobramycin")
Tobra16pairwise.m$Concentration <- c("16")

Tobra32 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "32")
Tobra32pairwise<-pairwise.t.test(Tobra32$OD, g=Tobra32$pH, p.adj="none")
Tobra32pairwise<-data.frame(Tobra32pairwise$p.value)
Tobra32pairwise$variable <- rownames(Tobra32pairwise)
Tobra32pairwise.m <- melt(Tobra32pairwise)
Tobra32pairwise.m <- na.omit(Tobra32pairwise.m)
Tobra32pairwise.m$Antibiotic <- c("Tobramycin")
Tobra32pairwise.m$Concentration <- c("32")

Tobra64 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "64")
Tobra64pairwise<-pairwise.t.test(Tobra64$OD, g=Tobra64$pH, p.adj="none")
Tobra64pairwise<-data.frame(Tobra64pairwise$p.value)
Tobra64pairwise$variable <- rownames(Tobra64pairwise)
Tobra64pairwise.m <- melt(Tobra64pairwise)
Tobra64pairwise.m <- na.omit(Tobra64pairwise.m)
Tobra64pairwise.m$Antibiotic <- c("Tobramycin")
Tobra64pairwise.m$Concentration <- c("64")


Mero0 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "0")
Mero0pairwise<-pairwise.t.test(Mero0$OD, g=Mero0$pH, p.adj="none")
Mero0pairwise<-data.frame(Mero0pairwise$p.value)
Mero0pairwise$variable <- rownames(Mero0pairwise)
Mero0pairwise.m <- melt(Mero0pairwise)
Mero0pairwise.m <- na.omit(Mero0pairwise.m)
Mero0pairwise.m$Antibiotic <- c("Meropenem")
Mero0pairwise.m$Concentration <- c("0")

Mero128 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "128")
Mero128pairwise<-pairwise.t.test(Mero128$OD, g=Mero128$pH, p.adj="none")
Mero128pairwise<-data.frame(Mero128pairwise$p.value)
Mero128pairwise$variable <- rownames(Mero128pairwise)
Mero128pairwise.m <- melt(Mero128pairwise)
Mero128pairwise.m <- na.omit(Mero128pairwise.m)
Mero128pairwise.m$Antibiotic <- c("Meropenem")
Mero128pairwise.m$Concentration <- c("128")

Mero256 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "256")
Mero256pairwise<-pairwise.t.test(Mero256$OD, g=Mero256$pH, p.adj="none")
Mero256pairwise<-data.frame(Mero256pairwise$p.value)
Mero256pairwise$variable <- rownames(Mero256pairwise)
Mero256pairwise.m <- melt(Mero256pairwise)
Mero256pairwise.m <- na.omit(Mero256pairwise.m)
Mero256pairwise.m$Antibiotic <- c("Meropenem")
Mero256pairwise.m$Concentration <- c("256")

Mero64 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "64")
Mero64pairwise<-pairwise.t.test(Mero64$OD, g=Mero64$pH, p.adj="none")
Mero64pairwise<-data.frame(Mero64pairwise$p.value)
Mero64pairwise$variable <- rownames(Mero64pairwise)
Mero64pairwise.m <- melt(Mero64pairwise)
Mero64pairwise.m <- na.omit(Mero64pairwise.m)
Mero64pairwise.m$Antibiotic <- c("Meropenem")
Mero64pairwise.m$Concentration <- c("64")

all.pairwise <- rbind(Gent0pairwise.m, Gent200pairwise.m, Gent400pairwise.m, Gent800pairwise.m,
                      Tobra0pairwise.m, Tobra16pairwise.m, Tobra32pairwise.m, Tobra64pairwise.m,
                      Mero0pairwise.m, Mero256pairwise.m, Mero64pairwise.m, Mero128pairwise.m)
all.pairwise$Strain <- STRAIN.name
all.pairwise.final <- rbind(all.pairwise, all.pairwise.final)

}

all.pairwise.final$adj <- p.adjust(all.pairwise.final$value, "bonf")
all.pairwise.final.filter<- subset(all.pairwise.final,adj < 0.05)
all.pairwise.final.filter$both <- paste(all.pairwise.final.filter$variable,all.pairwise.final.filter$variable.1, sep="_" )
all.pairwise.final.filter <- subset(all.pairwise.final.filter, !both=="basic_acidic")
'''

'''
all.fc.ttest.final <- c()
all.fc.ttest <- c()
temp<-c()
for (i in STRAIN.LIST) {
  STRAIN.name <- i
  temp <- subset(Antibiotics.FC, pH == "basic" | pH == "neutral")

  Gent200 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "200")
  Gent200pairwise<-pairwise.t.test(Gent200$FC,g=Gent200$pH,p.adj="none")
  Gent200pairwise<-data.frame(Gent200pairwise$p.value)
  Gent200pairwise$variable <- rownames(Gent200pairwise)
  Gent200pairwise.m <- melt(Gent200pairwise)
  Gent200pairwise.m <- na.omit(Gent200pairwise.m)
  Gent200pairwise.m$FC <- (t.test(Gent200$FC,p.adj="none", mu=1))["estimate"]
  Gent200pairwise.m$Antibiotic <- c("Gentamicin")
  Gent200pairwise.m$Concentration <- c("200")
  
  
  Gent400 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "400")
  Gent400pairwise<-pairwise.t.test(Gent400$FC,g=Gent400$pH,p.adj="none")
  Gent400pairwise<-data.frame(Gent400pairwise$p.value)
  Gent400pairwise$variable <- rownames(Gent400pairwise)
  Gent400pairwise.m <- melt(Gent400pairwise)
  Gent400pairwise.m <- na.omit(Gent400pairwise.m)
  Gent400pairwise.m$FC <- (t.test(Gent400$FC,p.adj="none", mu=1))["estimate"]
  Gent400pairwise.m$Antibiotic <- c("Gentamicin")
  Gent400pairwise.m$Concentration <- c("400")
  
  Gent800 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Gentamicin" & Concentration == "800")
  Gent800pairwise<-pairwise.t.test(Gent800$FC,g=Gent800$pH,p.adj="none")
  Gent800pairwise<-data.frame(Gent800pairwise$p.value)
  Gent800pairwise$variable <- rownames(Gent800pairwise)
  Gent800pairwise.m <- melt(Gent800pairwise)
  Gent800pairwise.m <- na.omit(Gent800pairwise.m)
  Gent800pairwise.m$FC <- (t.test(Gent800$FC,p.adj="none", mu=1))["estimate"]
  Gent800pairwise.m$Antibiotic <- c("Gentamicin")
  Gent800pairwise.m$Concentration <- c("800")
  
  Tobra16 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "16")
  Tobra16pairwise<-pairwise.t.test(Tobra16$FC,g=Tobra16$pH,p.adj="none")
  Tobra16pairwise<-data.frame(Tobra16pairwise$p.value)
  Tobra16pairwise$variable <- rownames(Tobra16pairwise)
  Tobra16pairwise.m <- melt(Tobra16pairwise)
  Tobra16pairwise.m <- na.omit(Tobra16pairwise.m)
  Tobra16pairwise.m$FC <- (t.test(Tobra16$FC,p.adj="none", mu=1))["estimate"]
  Tobra16pairwise.m$Antibiotic <- c("Tobramycin")
  Tobra16pairwise.m$Concentration <- c("16")

  Tobra32 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "32")
  Tobra32pairwise<-pairwise.t.test(Tobra32$FC,g=Tobra32$pH,p.adj="none")
  Tobra32pairwise<-data.frame(Tobra32pairwise$p.value)
  Tobra32pairwise$variable <- rownames(Tobra32pairwise)
  Tobra32pairwise.m <- melt(Tobra32pairwise)
  Tobra32pairwise.m <- na.omit(Tobra32pairwise.m)
  Tobra32pairwise.m$FC <- (t.test(Tobra32$FC,p.adj="none", mu=1))["estimate"]
  Tobra32pairwise.m$Antibiotic <- c("Tobramycin")
  Tobra32pairwise.m$Concentration <- c("32")
  
  Tobra64 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Tobramycin" & Concentration == "64")
  Tobra64pairwise<-pairwise.t.test(Tobra64$FC,g=Tobra64$pH,p.adj="none")
  Tobra64pairwise<-data.frame(Tobra64pairwise$p.value)
  Tobra64pairwise$variable <- rownames(Tobra64pairwise)
  Tobra64pairwise.m <- melt(Tobra64pairwise)
  Tobra64pairwise.m <- na.omit(Tobra64pairwise.m)
  Tobra64pairwise.m$FC <- (t.test(Tobra64$FC,p.adj="none", mu=1))["estimate"]
  Tobra64pairwise.m$Antibiotic <- c("Tobramycin")
  Tobra64pairwise.m$Concentration <- c("64")
  
  Mero64 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "64")
  Mero64pairwise<-pairwise.t.test(Mero64$FC,g=Mero64$pH,p.adj="none")
  Mero64pairwise<-data.frame(Mero64pairwise$p.value)
  Mero64pairwise$variable <- rownames(Mero64pairwise)
  Mero64pairwise.m <- melt(Mero64pairwise)
  Mero64pairwise.m <- na.omit(Mero64pairwise.m)
  Mero64pairwise.m$FC <- (t.test(Mero64$FC,p.adj="none", mu=1))["estimate"]
  Mero64pairwise.m$Antibiotic <- c("Meropenem")
  Mero64pairwise.m$Concentration <- c("64")
  
  Mero128 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "128")
  Mero128pairwise<-pairwise.t.test(Mero128$FC,g=Mero128$pH,p.adj="none")
  Mero128pairwise<-data.frame(Mero128pairwise$p.value)
  Mero128pairwise$variable <- rownames(Mero128pairwise)
  Mero128pairwise.m <- melt(Mero128pairwise)
  Mero128pairwise.m <- na.omit(Mero128pairwise.m)
  Mero128pairwise.m$FC <- (t.test(Mero128$FC,p.adj="none", mu=1))["estimate"]
  Mero128pairwise.m$Antibiotic <- c("Meropenem")
  Mero128pairwise.m$Concentration <- c("128")
  
  Mero256 <- subset(temp, Strain == paste(STRAIN.name) & Antibiotic=="Meropenem" & Concentration == "256")
  Mero256pairwise<-pairwise.t.test(Mero256$FC,g=Mero256$pH,p.adj="none")
  Mero256pairwise<-data.frame(Mero256pairwise$p.value)
  Mero256pairwise$variable <- rownames(Mero256pairwise)
  Mero256pairwise.m <- melt(Mero256pairwise)
  Mero256pairwise.m <- na.omit(Mero256pairwise.m)
  Mero256pairwise.m$FC <- (t.test(Mero256$FC,p.adj="none", mu=1))["estimate"]
  Mero256pairwise.m$Antibiotic <- c("Meropenem")
  Mero256pairwise.m$Concentration <- c("256")
  
  all.fc.ttest <- rbind(Gent200pairwise.m, Gent400pairwise.m, Gent800pairwise.m,
                        Tobra16pairwise.m, Tobra32pairwise.m, Tobra64pairwise.m,
                        Mero256pairwise.m, Mero64pairwise.m, Mero128pairwise.m)
  all.fc.ttest$Strain <- STRAIN.name
  all.fc.ttest.final <- rbind(all.fc.ttest, all.fc.ttest.final)
  
}

basic.fc.ttest<-all.fc.ttest.final
basic.fc.ttest$pH <- "basic"

acidic.fc.ttest<-all.fc.ttest.final
acidic.fc.ttest$pH <- "acidic"

all.fc.ttest.final <- rbind(acidic.fc.ttest, basic.fc.ttest)
all.fc.ttest.final$adj<- p.adjust(  all.fc.ttest.final$value, "bonf")
all.fc.ttest.final.filter<- subset(all.fc.ttest.final,adj < 0.05)

write.xlsx(all.fc.ttest.final.filter, 
           "/Volumes/GoogleDrive/My Drive/Stenotrophomonas/GrowthData/antibiotics.fc.sig.121118.xlsx")
'''

data.pH <- data[,-(8:10)]
data.pH.averages <- data.pH %>%
  group_by(Condition, Timepoint) %>%
  summarise(mean=mean(as.numeric(as.character(pH))))

pHplot<-ggplot()
pHplot <- pHplot + geom_point(data=data.pH,aes(x=as.numeric(as.character(Timepoint)), y=as.numeric(as.character(pH)), fill=Condition,
                                       shape=Condition),color="black",size=3)+
  scale_fill_manual(values=COLORS, name="")+
  scale_shape_manual(values=SHAPES, name="")+
  ylab("pH")+
  xlab("Hours")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=20),axis.text=element_text(size=16,colour="black"),
        legend.text = element_text(size=16),
        legend.title=element_text(size=20),
        legend.key = element_rect(colour = "white", fill="white"))+
  geom_line(data=data.pH.averages, size=1.5,aes(x=as.numeric(as.character(Timepoint)),y=mean, group=Condition, color=Condition 
  ))+
  scale_color_manual(values=COLORS,name="")


### Antibiotics growth curve
COLORS2 <-  c(Acidic="#e66101", Neutral="#cccccc", Basic="#5e3c99",
             unbuffered="black")


SHAPES2<-c(Acidic=22, Neutral=24, Basic=21,
          unbuffered=23)


require(grid) 
require(gridExtra)

tableA<-textGrob(c("a"), gp=gpar(fontsize=40), vjust=(-9))
tableB<-textGrob(c("b"), gp=gpar(fontsize=40), vjust=(0.3), hjust=(-1))
tableC<-textGrob(c("c"), gp=gpar(fontsize=40), vjust=(-9))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(ODplot)

 
PLOTOD <-arrangeGrob(arrangeGrob(tableA, (ODplot+theme(legend.position="none")), ncol=2, widths=c(1,12)),
             arrangeGrob(arrangeGrob(tableB, mylegend, nrow=2, heights=c(1,4)), antbxplot, ncol=2, widths=c(3,12)), 
             arrangeGrob(arrangeGrob(tableC,  FCplot, ncol=2, widths=c(1,12))),
             ncol=3, widths=c(1.5,3.75,3))

grid.arrange(PLOTOD)


ggsave("/Volumes/GoogleDrive/My Drive/Tara-KatrineGateway/Steno/StenoRevision2_12-11-18/FinalFigures/Fig2.eps", 
       device="eps", PLOTOD, width=18, height=9.3)

embed_fonts(file="./Fig2.eps",
            format="eps2write", outfile = "./Fig2.eps")

dev.off()


grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "blue", fill = NA))

acidsdata <- read.xlsx("/Volumes/GoogleDrive/My Drive/Stenotrophomonas/data/SmGrowthCurve_Acids.xlsx",1)
acidsdata <- na.omit(acidsdata)

COLORS2 <-  c(="#e66101", neutral="#cccccc", basic="#5e3c99",
              unbuffered="black" )


SHAPES2<-c(Acidic=22, Neutral=24, Basic=21,
           unbuffered=23)


acids.averages <- acidsdata %>%
  group_by(Timepoint, Condition) %>%
  summarise(mean=mean(as.numeric(as.character(OD))))

acidsplot <- ggplot()
acidsplot <- acidsplot+ geom_point(data=acidsdata, aes(x=Timepoint, y=OD, fill=Condition,shape=Condition),size=3)+
  #facet_wrap("Antibiotic", scales="free_x", ncol=1)+
 scale_fill_manual(values=c("#2E8B57", "#556B2F",
                            "#e66101", "#cccccc", "#5e3c99",
                            "#00FA9A", "black"), name="")+
 scale_shape_manual(values=c(25, 25, 22, 24, 21, 25,18), name="")+
  ylab("OD 600")+
  xlab("Incubation (h)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text=element_text(size=20),
        axis.title = element_text(size=20),axis.text=element_text(size=16,colour="black"),
        legend.text = element_text(size=16),
        legend.title=element_text(size=20),
        legend.key = element_rect(colour = "white", fill="white"))+
  geom_line(data=acids.averages, size=1.5,aes(x=Timepoint,y=mean, group=Condition, color=Condition
  ))+
  scale_color_manual(values=c("#2E8B57", "#556B2F",
                             "#e66101", "#cccccc", "#5e3c99",
                             "#00FA9A", "black"), name="")


+
 # scale_color_manual(values=COLORS, name="")

