##--------------------------------------------------------------------------------------------------#
## Variable correlation plot
##--------------------------------------------------------------------------------------------------#

library(corrplot)

#-#-# Correlation matrix #-#-#
filedir <- "E:/ProcessedData"
climData <- read.csv(paste0(filedir, "/ClimateData/bioclim_EWEMBI_1995_landonly.csv.gz"))
climData <- climData[,c("bio1","bio4","bio5","bio6","bio10","bio11","bio12","bio15","bio18","bio19")]
head(climData)

CorD <- cor(climData)
CorForP <- abs(CorD)

png("figures/correlationmatrix.png", width = 8, height = 7, units="in", res = 600)
corrplot(CorD, method="circle", bg = "white",addgrid.col = "gray10", 
         tl.col = "black",tl.cex = 0.8, p.mat = CorForP, sig.level = 0.7)
dev.off()
