##########################################################################
#               Future projections for GAM and GBM models                #
#                        LatLon 0.5 degree - global                      #
#              Predictions limited to nearest neighbour realm            #
#                             September 2017                             #
##########################################################################

# Load libraries
rm(list=ls())
packages <- c("PresenceAbsence", "plyr", "mgcv", "snowfall", "randomForest",
              "dplyr", "gbm", "matrixStats", "dismo")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

# Load GitHub package
#devtools::install_github("RS-eco/BioScen1.5_1")
#library(BioScen1.5_1)

# Set file directory
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir2 <- "/bigdata/mbiber" # shinichi
#filedir <- "/bigdata_local/mbiber" # ceremony - Mammals
#filedir <- "/home/mbiber/data" # ceremony - Birds
#filedir <- "/bigdata/mbiber/data" # darkstar
filedir <- "/media/matt/BS1p5_Amphibian"
#filedir <- "F:/"
filedir2 <- filedir

# Set working directory
#setwd("/scratch/home/mbiber/GitHub/BioScen1.5_1") #shinichi
#setwd("/home/mbiber/BioScen1.5_1") #ceremony
#setwd("/scratch/mbiber/BioScen1.5_1") #darkstar
setwd("/home/matt/Documents/Github/BioScen1.5_SDM")

########################################

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 1

# Model type
model_type <- c("GAM", "GBM", "MaxEnt", "RF")[2]

#-#-# Load in baseline and realm data #-#-#

## Load climate predictor combs
climCombs <- list(c("bio4", "bio5", "bio18", "bio19"),
                  c("bio4","bio5","bio12","bio15"),
                  c("bio4","bio5","bio12","bio15")) 
climCombs <- climCombs[i]

## Load climate data and blocks
climVar <- read.csv(
  paste0("data/Blocking_SU_1995_", 
         paste0(unlist(climCombs), collapse="_"), 
         ".csv"))[,c("x","y", unlist(climCombs), "block")]
colnames(climVar) <- c("x","y", unlist(climCombs),"block")

## Realm data 
#data(realm_coordinates) # Read in datafile with the realm coordinates
realm_coordinates <- read.csv("data/realm_coordinates.csv")

#-#-# Set filepath for results and future climate #-#-#

## Output file locations
timesteps <- c(1845, 1990, 1995, 2009, 2010, 2020, 2026, 2032, 2048, 2050, 2052, 2056, 
               2080, 2100, 2150, 2200, 2250)

# Save future projection files here
predPaths <- lapply(timesteps, function(x) 
  paste0(filedir, "/", taxa[i], "_", model_type, "_predictions_", x, "/"))

# Creates a folder if not already there
lapply(predPaths, function(x) if(!dir.exists(x)){dir.create(x)})

# Future climate data is here
future.data.root <- paste0(filedir, "/ClimateData/")
if(!dir.exists(future.data.root)){print("Climate data is missing!")}

# Original distribution of the species
bufferpath <- paste0(filedir, "/", taxa[i], "_Pseudoabsences/")

#-#-# Load baseline models and project #-#-#
## Model location
modelObjectLocation <- paste0(filedir2, "/", taxa[i], "_", model_type, "_Output/")

## Get all climate files (csv with xy and BioClimVar)
climAll <- list.files(paste0(future.data.root)) #All climate data
climAll <- climAll[!grepl(climAll, pattern="rcp85")] # Remove rcp85

#-#-# Get the realm neighbour grid #-#-#
Rneighbours <- read.csv("data/realm_neighbours.csv")
# Prepared csv file with neighbours for each realm (realm ID numbers)

#-#-# Set up the function to predict distributions #-#-# 
## Predict function
predict_func <- function(modelObjectLocation, model_type, spName, climData, clim){
  mod <- get(load(modelObjectLocation))
  lapply(mod, function(model){
    pred.block <- lapply(1:length(model),function(bk){
      if(model_type == "MaxEnt"){
        class(model[[bk]]$mod) <- "MaxEnt"
        PRED <- round(dismo::predict(model[[bk]]$mod, x=climData),4)
      }else if(model_type == "RF"){
        PRED <- round(predict(model[[bk]]$mod, newdata=climData, type="prob")[,2],4)
      } else{
        PRED <- round(predict(model[[bk]]$mod,newdata=climData, 
                              type="response",se.fit=FALSE),4)
      }
    })  
    names(pred.block) <- c(sapply(1:length(model),function(x){
      paste(clim,".",model[[x]]$block,sep="")
    }))
    return(pred.block)
  })
}

## Get list of models run and a list of future climate variables
modsAll <- list.files(modelObjectLocation, 
                      pattern = paste0("_model_output_", model_type),
                      full.names = T) # List all models

spNames <- sapply(basename(modsAll), function(x){
  paste0(strsplit(x, split="_")[[1]][1:2], collapse="_")
})
length(spNames)

AUC_data <- lapply(c("GAM", "GBM", "MaxEnt", "RF"), function(model_type){
  read.csv(paste0(filedir, "/AUCvalues_All_", 
                  model_type, "_", taxa[i], ".csv.xz"))})
AUC_data <- do.call(rbind, AUC_data)

# Aggregate the different AUC values from the 10 iterations per species
# and filter by AUC > 0.7
AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
  dplyr::summarise(mean = mean(AUC, na.rm=T)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
  group_by(Species, taxa) %>% dplyr::summarise(n = n()) %>% filter(n == 4)
spNames <- unique(spNames[spNames %in% AUC_sum$Species])
length(spNames)

## Existing files
modsMissing <- lapply(spNames, function(mods){
  if(any(!file.exists(sapply(predPaths, FUN=function(x){
    paste(x, "/", mods,"_PA_proj.csv.xz", sep="")})))){
    return(list.files(modelObjectLocation, pattern = mods, full.names = T))
  }
})
modsMissing <- Filter(Negate(is.null), modsMissing)
modsMissing <- unlist(modsMissing)
length(modsMissing); class(modsMissing)

# Check for corrupt files
#corrupt_files <- lapply(modsMissing, function(x){
#  data <- tryCatch(get(load(x)), error=function(e) e) #The prediction files
# if(inherits(data, "error")){
#    return(x)
#  } else{
#    return(NULL)
# }
#})
#corrupt_files <- unlist(Filter(Negate(is.null), corrupt_files))
#length(corrupt_files)
#file.remove(corrupt_files)

#-#-# Run prediction function for all species #-#-#
## Setup snowfall for parallisation
#(n <- ceiling(0.5*parallel::detectCores()))
#sfInit(parallel=TRUE, cpus=n)
#sfLibrary(mgcv); sfLibrary(gbm); sfLibrary(matrixStats); sfLibrary(randomForest);
#sfLibrary(dismo)#; sfLibrary(rJava)
#sfExport(list=c("predict_func", "climAll", "future.data.root", "modsAll", 
#                "modelObjectLocation", "timesteps", "model_type",
#                "climCombs", "predPaths", "bufferpath","realm_coordinates","Rneighbours"))
lapply(modsMissing[16:35], function(mods){
  spName <- paste(strsplit(basename(mods),split="_")[[1]][1:2],collapse="_") # Extract species name
  pseudoabsrep <- strsplit(basename(mods),split="_",fixed=TRUE)[[1]][3] # Get PA number
  if(any(!file.exists(sapply(predPaths, function(x) paste(x, "/", spName,"_",pseudoabsrep,"_proj.csv.xz", sep=""))))){ # Set to first time period folder
    possibleError <- tryCatch(getmodinfo <- get(load(mods)), error=function(e) e) # Get information on blocks and skip in case file is corupted
    if(!inherits(possibleError, "error")){
      ## What realm is the species in
      spbuffer <- get(load(paste0(bufferpath,spName,"_PA.Rdata")))[[1]] # Original distribution
      spbuffer <- subset(spbuffer,presence ==1) # Select the presence cells only
      spdist <- merge(spbuffer,realm_coordinates,all.x = TRUE) # Merge with realm data to get realm values
      RealmNr <- unique(spdist$Realm) # Unique realm values
      RealmNr <- RealmNr[!is.na(RealmNr)] # If extra cells in dist data, realm == NA thus additional cells get introduced
      
      if(length(RealmNr) > 0){
        Neighbours <- Rneighbours[which(Rneighbours$Realm %in% RealmNr),] # Make list of all neighbouring realms
        NeighValue <- as.list(Neighbours)
        NeighValue <- do.call(rbind,NeighValue)
        NeighValue <- unique(NeighValue[!is.na(NeighValue)])
        RealmPres <- realm_coordinates[(realm_coordinates$Realm %in% NeighValue),] # Get coordinates for all relevant realms
        RealmPres <- RealmPres[1:2]
      } else{
        # Some birds occur on big islands with no realm data!!!
        RealmPres <- spdist[1:2]
      }
      
      ## Run through all time periods
      lapply(1:length(timesteps), function(period){ # Climate files need to have projection year at the end (eg. 50, 80, 100) 
        if(!file.exists(paste0(predPaths[period], "/",spName,"_",pseudoabsrep,"_proj.csv.xz"))){
          climAllPeriod <- grep(climAll,pattern=timesteps[period],value=T) # Get all the climate files for the time period 
          
          ## Run through all GCM scenarios (rcps and models)
          spPredict <- lapply(climAllPeriod,function(clim){
            
            ## Get climate data
            climData <- read.csv(paste0(future.data.root,clim))[,c("x","y", unlist(climCombs))]
            climData <- merge(RealmPres, climData,all.x=TRUE)
            climName <- paste(strsplit(clim,split="_",fixed=TRUE)[[1]][2:3], collapse="_") # Depends on climate file names
            
            cat("\n","Starting: ", spName," ", model_type, " ", clim) # Print where I am at
            
            ## Run predictions
            predictData <- predict_func(modelObjectLocation=mods, 
                                        model_type=model_type, 
                                        spName=spName, 
                                        climData=climData, 
                                        clim=climName) # Predict function
            
            # Put predictions into dataframe
            outDATA <- as.data.frame(cbind(climData[,c("x","y")],predictData))
            
            # Define dataframe columns
            r <- paste0(model_type, "_",climName,"_block_",rep(1:(ncol(outDATA)-2)))
            colnames(outDATA) <- c("x", "y", r)
            
            # Remove data rows where all columns are 0
            outDATA <- outDATA[!!rowSums(abs(outDATA[-c(1,2)])),]
            return(outDATA)
          })
          all.mod <- Reduce(function(...) merge(...,by=c("x","y"),all.x=T),spPredict)
          readr::write_csv(all.mod, paste0(predPaths[period], "/",spName,"_",pseudoabsrep,"_proj.csv.xz"))
          gc()
          return(NULL)
        }
      })
    }
  }
})
sfStop()
q(save="no")
