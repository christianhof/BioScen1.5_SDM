# R Code for testing model performance with various variable combinations
# Written by Alke Voskamp & Matthias Biber

########################################

#' ## Install and load required packages

# Clear all memory
rm(list=ls())

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

# Specify file dir
#filedir <- "/bigdata/mbiber" # shinichi
#filedir <- "/bigdata_local/mbiber" # ceremony

# Set working directory
#setwd("/scratch/home/mbiber/GitHub/BioScen1.5_1") #shinichi
#setwd("/home/mbiber/BioScen1.5_1") # ceremony
filedir <- "M:/"
########################################

#' ## Read species data

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 3

########################################

#' ## Code to run GAMS for all species 
#' 
#' Including model calibration, distance sampling of absences and blocking
#' Modified to run for different variable combinations
#' Adapted from Naiara, David, Robbie and Christine 
#' August 2014

#' ## Load climate predictor combs and prepaire climate and block data 

# Different variables combinations
climCombs_3v <- as.list(as.data.frame(t(read.csv("data/VariableCombinations3_8.csv", header=T))))
climCombs_4v <- as.list(as.data.frame(t(read.csv("data/VariableCombinations4_8.csv", header=T))))
#climCombs_5v <- as.list(as.data.frame(t(read.csv("data/VariableCombinations5_8.csv", header=T))))
#climCombs <- c(climCombs_3v, climCombs_4v, climCombs_5v); rm(climCombs_3v, climCombs_4v, climCombs_5v)                                      
climCombs <- climCombs_5v

# Climate data
baseline <- read.csv("data/Blocking_SU_1995_bio1_bio4_bio12_bio15.csv")[-1]
head(baseline)

## Remove all unneeded columns
baseline <- baseline[,c("x","y",tolower(as.character(unique(unlist(climCombs)))),"block")]
head(baseline)

## Extract coordinates
coords <- baseline[,c("x","y")] #Get coordinates from climate data

#' ## Set the GAM model function that will be used to run the models later

#' GAM is different to the other models, 
#' it has internal cross validation (that's why it's so fast)

#' ## Set the modelling function up
source("R/GAM_func.R")

#' ## Run the GAMs including blocking and absence selection

# Specify input and output directory

sourceObs <- paste0(filedir, "/", taxa[i], "_Pseudoabsences/")
resultsPath <- paste0(filedir, "/", taxa[i], "_VariableSelectionModels_5v/")
if(!dir.exists(resultsPath)){dir.create(resultsPath)}

# Get all species pseudoabsence files
spFiles <- list.files(sourceObs, pattern=".Rdata", full.names=TRUE)

# Extract species names 
spNames <- lapply(spFiles,function(sp){
  name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse="_")})
})
spNames <- unlist(spNames)

# Load random subset of names
species <- read.csv(list.files("data/", pattern=paste0("random_sample_", tolower(taxa[i])), full.names=TRUE))

# Create subset list
spList <- spFiles[spNames %in% species$spName]

# Check number of missing Files
spMissing <- lapply(spList, function(sp){
  spname <- basename(sp)
  clim.var <- as.character(unlist(climCombs))
  climVarName <- paste(unlist(climCombs),collapse="_")
  sp <- strsplit(spname,split=".",fixed=T)[[1]][1]
  if(!file.exists(paste(resultsPath, sp,"_", paste(clim.var,collapse="_"), 
                        "_model_output_GAM.RData",sep=""))){
    return(sp)
  }
})
spMissing <- Filter(Negate(is.null), spMissing)
length(spMissing)

#Turn warning into error - if the model does not convert the code should stop rather 
#than giving a warning
options(warn=2)

# Set up snowfall to run parallel
sfInit(parallel=TRUE, cpus=ceiling(0.75*parallel::detectCores()))
sfLibrary(PresenceAbsence);sfLibrary(mgcv)
sfExport(list=c("GAM_eco","sourceObs","resultsPath","baseline","climCombs","spList")) 
#Import all the data, file path and model function needed to each CPU

# Run code for list
sfLapply(spMissing,function(sp){
  
  spname <- basename(sp)
  print(spname)
  
  lapply(climCombs, function(n){ #Loop through the different climate combinations
    
    clim.var <- tolower(as.character(unlist(n)))
    
    if(!file.exists(paste(resultsPath, "/", sp,"_", paste0(clim.var, collapse="_"), 
                          "_model_output_GAM_Eco_block.RData",sep=""))){ #Check if file exists already (to pick up where u stopped if code was )
      
      # Run model for each PA set
      mod <- lapply(1:10, function(y){
        PA <- paste0(spname, y) 
        
        spdata <- get(load(paste0(sourceObs,"/", spname, ".Rdata")))
        species.data <- spdata[[y]][,c("x","y","presence")]
        
        ## Select pseudo absence rep
        spPseudoRep <- na.omit(species.data) #Leave NAs out if there are any
        spPseudoRep <- merge(spPseudoRep,baseline,by=c("x","y"),all.x=T) #Merge the species data with the climate and blocking data
        spPseudoRep <- spPseudoRep[,c("x","y","presence","block",clim.var)] #Select only the relevant climate variables 
        spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep)) #Add ID column
        
        ## Summary stats
        ncells.pres <- nrow(spPseudoRep[spPseudoRep[,"presence"]==1,]) #Count the presences 
        
        if(ncells.pres >= 5){ #Skip restricted range species - this line is only to skip species below presence threshold (e.g. model on smaller grid)
          
          block.sum.all <- aggregate(presence~block,data=spPseudoRep,FUN=function(x)length(x)) #Sum the species data points per block
          block.sum.pres <- aggregate(presence~block,data=spPseudoRep,FUN=sum) #Sum the presences per block
          block.sum <- merge(block.sum.all,block.sum.pres,by="block",all.x=T)
          block.not.zero <- block.sum[block.sum[,3] > 0 & block.sum[,2] >= 10,] 
          num.block.not.zero <- nrow(block.not.zero) #Number of blocks that contain species presences
          
          ## Remove blocks that have zero presences
          if(as.numeric(ncells.pres) >= 5 & num.block.not.zero > 1){  # It has to be > 1 - There must be at least one block with presences and 5 presences overall
            
            block.include <- block.not.zero[,1]
            
            ## Model function from GAM
            GAM_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, PA=PA, 
                    clim.var=clim.var, fx=FALSE, k=-1, bs="tp",blocks=block.include)
          }
        } else{
          mod <- NULL
        }
      })
      mod <- Filter(Negate(is.null), mod)
      save(mod, file=paste(resultsPath, "/", sp,"_", paste(clim.var,collapse="_"), 
                           "_model_output_GAM.RData",sep=""), compress="xz")
    }
  })  
})
sfStop()

########################################

#' ## Get AUC stats across species and predictor combinations
#' July 2017

#rm(list=ls())

# List all files
GAMfiles <- list.files(paste0(filedir, "/", taxa[i], "_VariableSelectionModels/"))
head(GAMfiles)

spList <- lapply(GAMfiles,function(sp){
  name <- lapply(sp, function(x){paste(strsplit(x,"_",fixed=TRUE)[[1]][1:2],collapse="_")})
})

spList <- unlist(unique(spList))

# Try on one output
GAM1File <- get(load(paste0(filedir, "/", taxa[i], "_VariableSelectionModels/", GAMfiles[1])))
head(GAM1File)

## Look at AUC in first output 
AUC1 <- GAM1File[[1]]$AUC

# Set path to model files
mod.path <- paste0(filedir, "/", taxa[i], "_VariableSelectionModels/")
resultsPath <- paste0(filedir, "/", taxa[i], "_SummarisedModelOutput/")
if(!dir.exists(resultsPath)){dir.create(resultsPath)}

#modList <- GAMfiles[10]
#AUClist <- seq(1, 10, 1) 

#' Loop through all model output files to extract AUC and save summarized output
sfInit(parallel=TRUE, cpus=ceiling(0.55*parallel::detectCores()))
sfExport(list=c("GAMfiles","mod.path","resultsPath")) #Import all the data, file path and model function needed to each CPU

#Turn warning into error - if the model does not convert the code should stop rather than giving a warning
options(warn=2)

sfLapply(GAMfiles,function(modList){
  
  print(modList)
  
  ##Get species name 
  iteration <- strsplit(modList,split="_",fixed=TRUE)[[1]][3]
  
  # Check if model has 3 or 4 variables, by checking the length of the object
  LT <- length(unlist(strsplit(modList,split="_")))
  
  if(LT == 9){variables <- paste0(strsplit(modList,split="_",fixed=TRUE)[[1]][4:6],collapse="_")
  }else{variables <- paste0(strsplit(modList,split="_",fixed=TRUE)[[1]][4:7],collapse="_")}
  
  spName <- paste0(strsplit(modList,split="_",fixed=TRUE)[[1]][1:2],collapse="_")
  
  if(!file.exists(paste0(resultsPath,spName,"_",iteration,"_",variables,"_model_output_GAM.RData",sep=""))){
    
    ##Import species data
    mtry <- try(data.list <- get(load(paste0(mod.path,modList))))
    
    if(!inherits(mtry, "try-error")) {
      
      ##Extract AUC values
      AUCdata <- c()  
      
      for(i in 1:length(data.list)){
        blk <- data.list[i]
        AUC <- blk[[1]]$AUC
        print(AUC)
        AUCdata <- rbind(AUCdata,AUC)
      }
      
      finaldata <- c(spName, variables, iteration, unname(unlist(AUCdata)))
      finaldata <- as.data.frame(as.matrix(t(finaldata)))
      save(finaldata, file=paste(resultsPath,spName,"_",iteration,"_",variables,"_model_output_GAM.RData",sep=""),compress="xz") 
    }
    rm(data.list,AUCdata)
  }
}); sfStop()

# Read all results back in and save in one dataframe
All.files <- list.files(paste0(filedir, "/", taxa[i], "_SummarisedModelOutput"), 
                        pattern=".RData", full.names = TRUE)

all.sp.auc <- lapply(All.files,function(x){
  print(x)
  spdata <- get(load(x))
  head(spdata)
  return(spdata)
})
all.sp.auc <- do.call(rbind.fill,all.sp.auc)
colnames(all.sp.auc)<-c("Species","Variables","Iteration","AUC.1","AUC.2","AUC.3","AUC.4","AUC.5","AUC.6","AUC.7","AUC.8","AUC.9","AUC.10")
head(all.sp.auc)
nrow(all.sp.auc)
## Save AUC file as csv
write_csv(all.sp.auc, path=paste0(filedir, "/AUCvaluesAllModelsGAM_", taxa[i], ".csv")) 

# Aggregate the different AUC values from the 10 iterations per species
All.AUC <- read.csv(paste0(filedir, "/AUCvaluesAllModelsGAM_", taxa[i], ".csv"))
head(All.AUC)

AUCdata <- All.AUC[c(4:13)]
AUCdata <- rowMeans(AUCdata)
All.AUC$MeanAUCblocks <- AUCdata
All.Sub <- All.AUC[c(1,2,14)]

All.Sub.M <- tidyr::unite(All.Sub, newcol, c(Species, Variables),remove=FALSE)
All.Sub.M <- All.Sub.M[c(1,4)]
colnames(All.Sub.M) <- c("Models","AUC")
head(All.Sub.M)
nrow(All.Sub.M)

All.AUC.MeanPerSp <- aggregate(.~Models, data=All.Sub.M, mean)
head(All.AUC.MeanPerSp)
nrow(All.AUC.MeanPerSp)

readr::write_csv(All.AUC.MeanPerSp, path=paste0(filedir,"/AUCvaluesAllModelsGAM_Sum_", taxa[i], ".csv"))

# Rank combinations from aggregated file
All.AUC <- read.csv(paste0(filedir, "/AUCvaluesAllModelsGAM_Sum_", taxa[i], ".csv"))
head(All.AUC)

# Spit Models
All.AUC$Species <- sapply(All.AUC$Models, function(x) paste0(strsplit(as.character(x), split="_", fixed=TRUE)[[1]][1:2], collapse="_")) 
All.AUC$Models <- sapply(All.AUC$Models, function(x) paste0(strsplit(as.character(x), split="_", fixed=TRUE)[[1]][-c(1,2)], collapse="_")) 

##Loop through species and rank by AUC value
spList <- unique(as.vector(All.AUC$Species))

AUCdata <- lapply(spList,function(x){
  print(x)
  sp <- subset(All.AUC,Species == x)
  sp <- sp[with(sp,order(-AUC)),]
  rankNr <- c(1:nrow(sp))
  sp$rank <- rankNr
  return(sp)
})

FinalRank <- do.call(rbind,AUCdata)
head(FinalRank)

readr::write_csv(FinalRank, path=paste0(filedir, "/FinalRank_", taxa[i], ".csv"))
