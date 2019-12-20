#-#-# Creat csv with all possible climate combinations three and four variables #-#-#

rm(list=ls(all=TRUE))
library(gtools)
filedir <- "E:/ProcessedData"
climvar <- read.csv(paste0(filedir, "/ClimateData/bioclim_EWEMBI_1995_landonly.csv.gz"))
head(climvar)
length(climvar)
AllVariables<-c("bio1","bio4","bio5","bio6","bio10","bio11","bio12","bio15","bio18","bio19")
VariableCombinations <- expand.grid(AllVariables,AllVariables,AllVariables)
head(VariableCombinations)
nrow(VariableCombinations)
VariableCombinations$uniN <- apply(VariableCombinations,1,function(x) length(unique(c(x)))) #how many unique terms are in the combination
VariableCombinations <- subset(VariableCombinations, uniN == 3)[,1:3]
nrow(VariableCombinations)
VariableCombinations <- as.data.frame(t(apply(VariableCombinations,1,function(x) mixedsort(x)))) #sort each row by number 
VariableCombinations <- VariableCombinations[!duplicated(VariableCombinations), ] #remove all duplicate rows
# Need to manually remove correlation combinations!!!
#write.csv(VariableCombinations, "data/VariableCombinations3_8.csv", row.names=F)

filedir <- "E:/ProcessedData"
climvar <- read.csv(paste0(filedir, "/ClimateData/bioclim_EWEMBI_1995_landonly.csv.gz"))
head(climvar)
length(climvar)
AllVariables<-c("bio1","bio4","bio5","bio6","bio10","bio11","bio12","bio15","bio18","bio19")
VariableCombinations <- expand.grid(AllVariables,AllVariables,AllVariables,AllVariables)
head(VariableCombinations)
nrow(VariableCombinations)
VariableCombinations$uniN <- apply(VariableCombinations,1,function(x) length(unique(c(x)))) #how many unique terms are in the combination
VariableCombinations <- subset(VariableCombinations, uniN == 4)[,1:4]
nrow(VariableCombinations)
VariableCombinations <- as.data.frame(t(apply(VariableCombinations,1,function(x) mixedsort(x)))) #sort each row by number 
VariableCombinations <- VariableCombinations[!duplicated(VariableCombinations), ] #remove all duplicate rows
# Need to manually remove correlation combinations!!!
#write.csv(VariableCombinations, "data/VariableCombinations4_8.csv", row.names=F)

filedir <- "E:/ProcessedData"
climvar <- read.csv(paste0(filedir, "/ClimateData/bioclim_EWEMBI_1995_landonly.csv.gz"))
head(climvar)
length(climvar)
AllVariables<-c("bio1","bio4","bio5","bio6","bio10","bio11","bio12","bio15","bio18","bio19")
VariableCombinations <- expand.grid(AllVariables,AllVariables,AllVariables,AllVariables,AllVariables)
head(VariableCombinations)
nrow(VariableCombinations)
VariableCombinations$uniN <- apply(VariableCombinations,1,function(x) length(unique(c(x)))) #how many unique terms are in the combination
VariableCombinations <- subset(VariableCombinations, uniN == 5)[,1:5]
nrow(VariableCombinations)
VariableCombinations <- as.data.frame(t(apply(VariableCombinations,1,function(x) mixedsort(x)))) #sort each row by number 
VariableCombinations <- VariableCombinations[!duplicated(VariableCombinations), ] #remove all duplicate rows
VariableCombinations <- as.data.frame(apply(VariableCombinations, 2, toupper))
# Need to manually remove correlation combinations!!!
write.csv(VariableCombinations, "data/VariableCombinations5_8.csv", row.names=F)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#