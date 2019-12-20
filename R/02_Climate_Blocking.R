#' Create model blocks according to climate variable
 
########################################

#' ## Install and load required packages

# Automatically install required packages, which are not yet in library
packages <- c("base", "lattice", "randomForest", "sp", "spatstat", "celestial", "maptools", "stats", 
              "graphics", "parallel", "utils", "mgcv", "deldir", "SDMTools", "raster", "rgdal",
              "snowfall", "ggplot2", "blockTools", "PresenceAbsence", 
              "zoo", "plyr", "dplyr", "tidyr", "readr")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

########################################

#' ## Set file directory

filedir <- "E:/ProcessedData" # Desktop

########################################

#' ## Read species data

# Set taxa
taxa <- c("Amphibian", "Ter_Mammal", "Ter_Bird")

#' Land coordinates
land <- read.csv("extdata/Realm_coordinates.csv")
land$Realm <- as.factor(land$Realm)
ggplot()+geom_raster(data=land,aes(x=x,y=y,fill=Realm))

########################################

#' ## Blocking the global climate data by ecoregion 
#' Based on Robies' blocking method GCB paper 2013 
#' July 2014

#' See create_blockingsampleunits.R for creating the blockingsampleunits.csv file.

#' ## Create the blocks using the baseline climate data

#' Read in the data

# Read in the csv with the Blocking Samples
sample.units.id<- read.csv("data/blockingsampleunits.csv")
colnames(sample.units.id) <- c("x","y","id.sample")

# Read in the climate data
filedir <- "E:/ProcessedData"
climData <- read.csv(paste0(filedir, "/ClimateData/bioclim_EWEMBI_1995_landonly.csv.gz"))
names(climData)
out<-climData
sample.units.id <- merge(sample.units.id,out,by=c("x","y"),all.x=FALSE)
names(sample.units.id)
head(sample.units.id)

#' Do the blocking

# !!Choose your climate variables you want to model with and block with those!!
var_combs <- list(c("bio1", "bio4", "bio12", "bio15"), 
                  c("bio4", "bio5", "bio12", "bio15"),
                  c("bio4", "bio5", "bio18", "bio19"))

# Aggregate to sample regions 
for(i in 1:length(var_combs)){
  sample.unit.mean <- aggregate(sample.units.id[,unlist(var_combs[i])], by=list(sample.units.id$id.sample), mean)
  colnames(sample.unit.mean) <- c("id.sample", sapply(var_combs[i], FUN=function(x) paste0(x, "m")))
  head(sample.unit.mean)
  sample.unit.var <- aggregate(sample.units.id[,unlist(var_combs[i])], by=list(sample.units.id$id.sample), var)
  colnames(sample.unit.var) <- c("id.sample", sapply(var_combs[i], FUN=function(x) paste0(x, "v")))
  head(sample.unit.var)
  sample.unit.var[is.na(sample.unit.var)] <- 0
  sample.unit.all <- merge(sample.unit.mean,sample.unit.var,by="id.sample",all.x=TRUE)
  length(sample.unit.all[,1])
  # names(sample.unit.all)
  head(sample.unit.all)
  
  #Create blocks, dividing polygons into orthogonal blocks on the basis of the climate data.
  blocks <- blockTools::block(sample.unit.all, n.tr=10, id.vars='id.sample') ## is this where you decide the number of blockds to use i.e. n.tr=5 or 10
  blocks <- blockTools::assignment(blocks, namesCol=as.character(1:10))$assg[[1]][1:10] ## assign subpols to one of 10 blocks 
  blocks <- Reduce(rbind, mapply(function(id.sample, block) data.frame(id.sample, block), id.sample=blocks, block=as.list(1:10),SIMPLIFY=F) )## turn this into a dataframe 
  blocks <- blocks[!is.na(blocks$id.sample),] ## remove any polygons that haven't been assigned to a block (we deal with this later).
  head(blocks)
  clidat <- merge(sample.unit.all, blocks, by=c('id.sample'), all=T) 
  head(clidat)
  clidat$block[is.na(clidat$block)] <- sample(1:10,1)
  clidat <- clidat[,c(1,10)] # id sample and block ## this value will need to change depending on how many variables I have
  names(clidat)
  
  t <- merge(sample.units.id, clidat,by="id.sample",all.x=TRUE)
  head(t)
  lattice::levelplot(block~x+y,data=t)
  names(sample.units.id)
  write.csv(t, paste0("data/Blocking_SU_1995_", paste(unlist(var_combs[i]), collapse="_"), ".csv"),
            row.names=F)
}
