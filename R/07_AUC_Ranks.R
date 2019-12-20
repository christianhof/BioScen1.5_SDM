#-#-# Stacked barchart 400 species data #-#-#
rm(list=ls())

library(data.table)
library(plyr)
library(ggplot2)
library(lattice)
library(fields)
library(maptools)
library(colorRamps)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(reshape2)

# Specify file dir
filedir <- "C:/ProcessedData" # Desktop

# Set taxa
taxa <- c("Amphibian", "Ter_Mammal", "Ter_Bird")

for(i in c(1,2,3)){
  #-#-# Summarize frequenzies and plot as stacked bar chart #-#-#
  AUCall <- read.csv(paste0(filedir, "/FinalRank_", taxa[i], ".csv"))
  
  head(AUCall)
  nrow(AUCall)
  
  #-#-# Select the ten best variable combinations to display in graph #-#-#
  AUCallTop <- AUCall ## Add data here
  AUCallTop$rank[AUCallTop$rank >= 4] <- 0
  
  head(AUCallTop)
  
  AUCallTop <- AUCallTop[,c("Species","Models","rank")]
  head(AUCallTop)
  
  AUCTopTable<-data.frame(table(AUCallTop$Models, AUCallTop$rank))
  colnames(AUCTopTable) <- c("ClimateVariable","Rank","Frequency")
  head(AUCTopTable)
  
  AUCTopTable$Rank <- as.numeric(as.character(AUCTopTable$Rank))
  AUCTopTable <- subset(AUCTopTable,Rank >= 1)
  head(AUCTopTable)
  #View(AUCTopTable)
  
  AUCTopTable$Frequency <- as.numeric(as.character(AUCTopTable$Frequency))
  AUCTopTable <- aggregate(Frequency ~ ClimateVariable, AUCTopTable, sum)
  AUCTopTable <- AUCTopTable[order(-AUCTopTable$Frequency),]
  head(AUCTopTable)
  
  AUCTopVariables <- AUCTopTable[1:10,]
  
  TopVariableList <- as.vector(AUCTopVariables$ClimateVariable)
  
  #-#-# Subset the entire results data frame choosing only the best variables #-#-#
  AUCall <- AUCall
  
  AUCall$rank[AUCall$rank >= 10] <- "Other"
  head(AUCall)
  nrow(AUCall)
  
  AUCSub <- AUCall[,c("Species","Models","rank")]
  head(AUCSub)
  
  AUCFreqTable<-data.frame(table(AUCSub$Models, AUCSub$rank))
  colnames(AUCFreqTable) <- c("ClimateVariable","Rank","Frequency")
  head(AUCFreqTable)
  
  AUCallFinal <- AUCFreqTable[AUCFreqTable$ClimateVariable %in% TopVariableList, ]
  library(dplyr)
  AUCallFinal %>% arrange(Rank, desc(Frequency)) %>% head()
  nrow(AUCallFinal)
  #View(AUCallFinal)
  AUCallFinal <- subset(AUCallFinal, Frequency > 0)
  
  #-#-# Set colour scheme #-#-#
  PaletteBlue2 <-c('blue4','dodgerblue4','deepskyblue','gray20','gray28','gray39','gray49','gray53','gray63','gray73')
  
  #-#-# Extract all label names #-#-#
  #labellist <- as.vector(AUCFreqTable$ClimateVariable)
  #labellist <- unique(labellist)
  
  testMax <- (nrow(AUCall))/23
  
  #-#-# Plot the variables #-#-#
  p1s <- ggplot(AUCallFinal, aes(x = ClimateVariable, y = Frequency)) +
    geom_bar(aes(fill = Rank), stat="identity") +
    scale_fill_manual(values=PaletteBlue2)+
    scale_y_continuous() +
    guides(fill = guide_legend(ncol = 1))+
    theme(panel.background=element_rect(fill="white",colour="white"),
          panel.grid=element_blank(),
          plot.title = element_text(lineheight=2, face="bold",hjust = 0),
          axis.text=element_text(size=8, colour="black"), 
          axis.title=element_text(size=8),
          axis.line=element_line(colour="black"))+
    theme(axis.text.x=element_text(angle = -90, hjust = 0))
  
  print(p1s)
  
  grobframe <- arrangeGrob(p1s, ncol = 1, nrow=1)
  plot(grobframe)
  
  tiff(file = paste0("figures/StackedBar_", taxa[i], "_Top10.tiff"), width = 9000, 
       height = 4000, units = "px",  compression = "lzw", res = 800) 
  plot(grobframe) 
  dev.off() 
  
  #grobframe <- arrangeGrob(p1,p2,p3, ncol = 3, nrow=1)
  
  #-#-#
  
  boxplot(AUCall$AUC~AUCall$Models,data=AUCall)
}
