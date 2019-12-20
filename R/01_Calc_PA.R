#' # Calculate the pseudo absence for each species

########################################

#' ## Install and load required packages

# Clear all memory
rm(list=ls())

# Automatically install required packages, which are not yet in library
packages <- c("base", "lattice", "sp", "spatstat", "maptools", "stats", 
              "SDMTools", "raster", "readr")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

########################################

#' ## Set file directory

# Specify file dir
filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/scratch/mbiber/data" # darkstar/ceremony
#filedir <- "E:/ProcessedData" # Desktop

########################################

#' Get species names
data(amphibians, package="rasterSp")
data(ter_birds, package="rasterSp")
data(ter_mammals, package="rasterSp")

# Set taxa
taxa <- c("Amphibian", "Ter_Mammal", "Ter_Bird")

# Read in the path to species files
spfilePathWGS <- paste0(filedir, "/SpeciesData/")

# Check which files are already there 
available_files <- list.files(spfilePathWGS)
available_names <- sapply(available_files, FUN=function(x) 
  strsplit(as.character(x), split=".tif")[[1]][1])
rm(available_files)

# Only use species for which raster files exist
speciesList <- list()
speciesList[[1]] <- amphibians$binomial[amphibians$binomial %in% available_names]
speciesList[[2]] <- ter_mammals$binomial[ter_mammals$binomial %in% available_names]
speciesList[[3]] <- ter_birds$SCINAME[ter_birds$SCINAME %in% available_names]

########################################

#' ## Calculate pseudo absences
#'
#' Calculation is based on the distance to a species' range 
#' using Naiara's distance weighting "One over distance 2"

# Run code for all three taxa
for(i in 1:length(taxa)){
  
  # Read in the path to result files
  resultspath <- paste0(filedir, "/", taxa[i], "_Distances/")
  
  if(!dir.exists(resultspath)){dir.create(resultspath)}
  
  # Load paths of raster files
  sp.path <- lapply(speciesList[[i]],function(x){
    species <- paste0(spfilePathWGS, x,".tif")
  })
  sp.path <-  do.call(rbind,sp.path)
  sp.path[[1]]
  
  # Initialise parallel processing
  sfInit(parallel=TRUE, cpus=detectCores()-1)
  sfLibrary(spatstat); sfLibrary(sp); sfLibrary(raster); 
  sfLibrary(maptools); sfLibrary(SDMTools)
  
  # Source the distance.calc function
  source("R/distance_func.R")
  
  # Import all the data and data paths needed to each CPU
  sfExport(list=c("resultspath", "sp.path", "distance.calc", 
                  "spfilePathWGS", "land")) 
  
  # Run distance.calc function parallel
  system.time(
    sfLapply(sp.path,function(x) distance.calc(x))
    #lapply(sp.path[1:5],function(x) distance.calc(x))
  )
  # Stop clusters
  sfStop()
  
  # Check output
  test <- get(load(list.files(resultspath, full.names=TRUE)[1]))
  head(test)
  ggplot() + geom_raster(data=test,aes(x=x,y=y,fill=OneOverDist2))
  
  ## Select pseudo absences based on the distance to a species' range
  
  # Set file path and get distance data
  spDistDir <- resultspath
  spName <- list.files(spDistDir)
  spPresDir <- spfilePathWGS # paste0(filedir, "/SpeciesData/")
  spName[[1]]
  
  # Specify output file dir
  filetest <-  paste0(filedir, "/", taxa[i], "_Pseudoabsences/")
  
  # Create file dir if necessary
  if(!dir.exists(filetest)){dir.create(filetest)}
  
  # Initialise parallel processing
  sfInit(parallel=TRUE, cpus=detectCores()-1)
  sfLibrary(base);sfLibrary(lattice);sfLibrary(raster);
  
  # Source the distance.calc function
  source("R/PA_func.R")
  
  # Import all the data and data paths needed to each CPU
  sfExport(list=c("spDistDir", "spName", "spPresDir", "filetest", "PA.calc")) 
  
  system.time(
    sfLapply(spName,function(sp) PA.calc(sp))
  )
  sfStop()
  
  # Check the output
  test <- get(load(list.files(filetest, full.names=TRUE)[1]))
  head(test[["PA1"]])
  ggplot()+geom_raster(data=test,aes(x=x,y=y,fill=factor(presence)))
}