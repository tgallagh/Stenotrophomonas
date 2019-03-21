#### smear plots
### 07/28/17

require(xlsx)
require(ggplot2)


all57 <- read.xlsx("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/all.5vs7.genes.xlsm", 1)
all57$Threshold <- factor(all57$Threshold, levels = c("Not sig", "FDR", "FDR and FC"))
plot <- ggplot()
plot + geom_point(data = all57, aes(x=logCPM, y=logFC, color = Threshold))+
  scale_color_manual(values = c("black", "black", "#e6550d"))+
  scale_shape_manual(values= c(17, 20))+    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(),        
legend.position ="none",
axis.text=element_text(size=20, color="black"),
axis.title=element_text(size=20),
axis.line = element_line(colour = "black"))+
  labs(y=expression(paste(log["2"],FC)))+
  labs(x=expression(paste(log["2"],CPM)))

all79 <- read.xlsx("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/all.7vs9.genes.xlsx", 1)

all79 <- read.delim("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/all.7vs9.genes.txt")

all79$Threshold <- factor(all79$Threshold, levels = c("Not sig", "FDR", "FDR and FC"))
plot<-ggplot()

plot <- ggplot()
plot + geom_point(data = all79, aes(x=logCPM, y=logFC, color = Threshold))+
  scale_color_manual(values = c("black", "black", "purple"))+
  scale_shape_manual(values= c(17, 20))+    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),        
      legend.position="none",
        axis.text=element_text(size=20, color="black"),
        axis.title=element_text(size=20),
        axis.line = element_line(colour = "black"))+
  labs(y=expression(paste(log["2"],FC)))+
  labs(x=expression(paste(log["2"],CPM)))

