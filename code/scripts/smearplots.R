#### smear plots
### 07/28/17

require(xlsx)
require(ggplot2)


all57 <- read.xlsx("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/all.5vs7.genes.xlsm", 1)
plot <- ggplot()
plot + geom_point(data = all57, aes(x=logCPM, y=-logFC, color = sig, shape = method))+
  scale_color_manual(values = c("#D31E1E", "#C46363", "black"))+
  scale_shape_manual(values= c(17, 20))+    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("logFC") +
  xlab("logCPM")+
  ggtitle("pH 5 differentially expressed genes")
  
all79 <- read.xlsx("~/GoogleDrive/Stenotrophomonas/data/processed/edgeR/all.7vs9.genes.xlsm", 1)
plot <- ggplot()
plot + geom_point(data = all79, aes(x=logCPM, y=-logFC, color = sig, shape = method))+
  scale_color_manual(values = c("#00078C", "#5E65DA", "black"))+
  scale_shape_manual(values= c(17, 20))+    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("logFC") +
  xlab("logCPM")+
  ggtitle("pH 9 differentially expressed genes")
