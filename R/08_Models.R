#############################################################################################
#                            Code to run models for all species                             #
#                                    LatLon data 0.5 degree                                 #
#               Blocking by Ecoregions if species have 50 or more presence                  #
#                          otherwise 30/70 split for model validation                       #
#                        Code adapted from David, Robbie and Naiara                         #
#############################################################################################

rm(list=ls())

# Load packages
# Automatically install required packages, which are not yet in library
packages <- c("mgcv", "PresenceAbsence", "snowfall", "ggplot2", "dplyr",
              "gbm", "dismo", "randomForest")#, "rJava")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

#' ## Set file directory

# Specify file dir
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/bigdata_local/mbiber" # ceremony - Mammals
#filedir <- "/home/mbiber/data" # Ceremony - Birds
filedir <- "/bigdata/mbiber/data" # darkstar

# Set working directory
#setwd("/scratch/home/mbiber/GitHub/BioScen1.5_1") # shinichi
#setwd("/home/mbiber/BioScen1.5_1") # ceremony
setwd("/scratch/mbiber/BioScen1.5_1") #darkstar

########################################

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 3

# Set model_type
model_type <- c("GAM", "GBM", "MaxEnt", "RF")[4]

########################################

#-#-# Load the baseline climate #-#-#
# Load climate predictor combs
climCombs <- list(c("bio4", "bio5", "bio18", "bio19"),
                  c("bio4","bio5","bio12","bio15"),
                  c("bio4", "bio5", "bio12", "bio15")) # Climate variables for the models make sure its the right combination
climCombs <- climCombs[i]

## Load climate data and blocks
climVar <- read.csv(paste0("data/Blocking_SU_1995_", 
                           paste0(unlist(climCombs), collapse="_"), ".csv"))[,c("x","y", unlist(climCombs), 
                                                                                "block")]
colnames(climVar) <- c("x","y", unlist(climCombs),"block")

########################################

#' ## Read species data

#-#-# Set the file paths #-#-#
sourceObs <- paste0(filedir, "/", taxa[i], "_Pseudoabsences")
resultsPath <- paste0(filedir, "/", taxa[i], "_", model_type, "_Output")
if(!dir.exists(resultsPath)){dir.create(resultsPath)}
if(model_type == "GBM"){
  plotPath <- paste0(filedir, "/", taxa[i], "_", model_type, "_plots")
  if(!dir.exists(plotPath)){dir.create(plotPath)}
} else{
  plotPath <- NA
}

#-#-# List the species files #-#-#
spFiles <- list.files(sourceObs, full.names=TRUE)
head(spFiles)

########################################

#Turn warning into error -
#if the model does not convert the code should stop rather than giving a warning
options(warn=2)

#-#-# Set the model function that will be used to run the models later #-#-#
source("R/GAM_func.R")
source("R/GBM_func.R")
source("R/MaxEnt_func.R")
source("R/RF_func.R")

# Read data with number of presences per species
#library(dplyr)
#sp_split <- read.csv("extdata/no_records_species_groups.csv")
#taxa_long <-  c("Amphibians", "Terrestrial Mammals", "Terrestrial Birds")
#sp_split <- sp_split %>% filter(group == taxa_long[i]) %>% 
#  filter(sum >= 50) %>% select(species)
#sp_split <- sp_split %>% filter(group == taxa_long[i]) %>%  
#filter(sum >= 10 & sum < 50) %>% select(species)

# Subset files according to correct number of presences
#spFiles <- spFiles[lapply(basename(spFiles), function(x) 
#  paste(strsplit(x, split="_")[[1]][1:2], collapse=" ")) %in% sp_split$species]

# Remove existing files from list
spMissing <- lapply(spFiles, function(sp){
  spname <- basename(sp)
  clim.var <- as.character(unlist(climCombs))
  climVarName <- paste(unlist(climCombs),collapse="_")
  sp <- strsplit(spname,split=".",fixed=T)[[1]][1]
  if(!file.exists(paste(resultsPath, "/", sp, "_",climVarName,
                        "_model_output_", model_type, "_Eco_block.RData",sep=""))){
    if(!file.exists(paste(resultsPath, "/", sp, "_",climVarName,
                          "_model_output_", model_type, "_30_70.RData",sep=""))){
      if(!file.exists(paste(resultsPath, "/", sp, "_",climVarName,
                            "_model_output_", model_type, "_30_70_MissingEco.RData",sep=""))){
        return(sp)
      }
    }
  }
})
spMissing <- Filter(Negate(is.null), spMissing)
length(spMissing)

# Read data
AUC_data <- lapply(c("GAM", "GBM", "MaxEnt", "RF"), function(model_type){
    read.csv(paste0(filedir, "/AUCvalues_All_", 
                    model_type, "_", taxa[i], ".csv.xz"))})
AUC_data <- do.call(rbind, AUC_data)

# Aggregate the different AUC values from the 10 iterations per species
# and filter by AUC > 0.7
AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
  summarise(mean = mean(AUC, na.rm=T)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
  group_by(Species, taxa) %>% summarise(n = n()) %>% filter(n == 4)

spNames <- sub("_PA", "", spMissing)
spMissing <- unique(spMissing[spNames %in% AUC_sum$Species])
length(spMissing)

# Set up snowfall to run the GAMs including blocking/30-70 split and absence selection #
library(snowfall)
sfInit(parallel=TRUE, cpus=ceiling(0.05*parallel::detectCores()))
sfLibrary(PresenceAbsence); sfLibrary(mgcv); sfLibrary(gbm); 
sfLibrary(dismo); sfLibrary(dplyr); sfLibrary(randomForest)
#sfLibrary(rJava)

# Import all the data and data paths needed to each CPU
sfExport(list=c("GAM_split", "GAM_eco","sourceObs", "resultsPath", 
                "climCombs", "climVar", "GBM_eco", "GBM_split", "model_type",
                "plotPath", "RF_eco", "RF_split", "MaxEnt_eco", "MaxEnt_split")) 

# Run code
source("R/model_run.R")
sfLapply(spMissing, model_run)
sfStop()
# system('shutdown -s')

# Test output
#sp <- "Spinomantis_peraccae_PA"
#mod <- get(load(paste(resultsPath, "/", sp, "_",paste(unlist(climCombs),collapse="_"),
#               "_model_output_", model_type, "_Eco_block.RData",sep="")))