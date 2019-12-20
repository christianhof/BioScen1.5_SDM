#-#-# Add dispersal to predictions data

# Load libraries
rm(list=ls())
packages <- c("raster", "readr", "rgeos", "tidyr", "sp", "ggplot2", "snowfall")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

########################################

# Set file directory
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/bigdata_local/mbiber" # ceremony - Mammals
#filedir <- "/home/mbiber/data" # ceremony - Birds
filedir <- "/bigdata/mbiber/data" # Darkstar

# Set working directory
#setwd("/scratch/home/mbiber/GitHub/BioScen1.5_1") #shinichi
#setwd("/home/mbiber/BioScen1.5_1") #ceremony
setwd("/scratch/mbiber/BioScen1.5_1") #darkstar

########################################

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 3

# Model type
k <- 4; model_type <- c("GAM", "GBM", "MaxEnt", "RF")[k]

#-#-# Get future data #-#-#
predpath <- paste0(filedir, "/", taxa[i], "_", model_type, "_predictions/")

# List all files
modsAll <- list.files(predpath, full.names = T)
spNames <- sapply(basename(modsAll), function(x){
  paste0(strsplit(x, split="_")[[1]][1:2], collapse="_")
})
length(spNames)

# Read AUC data
AUC_data <- lapply(c("GAM", "GBM", "MaxEnt", "RF"), function(model_type){
  readr::read_csv(paste0(filedir, "/AUCvalues_All_", 
                         model_type, "_", taxa[i], ".csv"))})
AUC_data <- do.call(rbind, AUC_data)

# Aggregate the different AUC values from the 10 iterations per species
# and filter by AUC > 0.7
library(dplyr)
AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
  summarise(mean = mean(AUC, na.rm=T)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
  group_by(Species, taxa) %>% summarise(n = n()) %>% filter(n == 4)

spNames <- unique(spNames[spNames %in% AUC_sum$Species])
length(spNames)

#-#-# Set filepath original distribution and future #-#-#
SpeciesData <- paste0(filedir, "/SpeciesData/")
resultspath <- paste0(filedir, "/", taxa[i], "_", model_type, "_results_climate/")
if(!dir.exists(resultspath)){dir.create(resultspath)}

#-#-# Get area of each gridcell #-#-#
area <- read.csv("data/realm_coordinates.csv")
area <- area[,c("x", "y", "area")]
colnames(area) <- c("x","y","areaKM2")

# Get missing names
names_mis <- lapply(spNames, function(x){
  if(!file.exists(paste0(resultspath,  x, "_", model_type, "_dispersal.csv.xz"))){
    return(x)
  }
})
names_mis <- Filter(Negate(is.null), names_mis)
length(names_mis)

# List all summary prediction files, where dispersal is missing
#files <- unlist(lapply(names_mis, function(species){
#  list.files(paste0(filedir, "/", taxa[i], "_", model_type, "_predictions"), 
#             pattern=species, full.names=T)}))
#length(files)

# Check for corrupt files
#corrupt_files <- lapply(files, function(x){
#  data <- tryCatch(readr::read_csv(x), error=function(e) e) #The prediction files
#  if(inherits(data, "error")){
#    return(x)
#  } else{
#    return(NULL)
#  }
#})
#corrupt_files <- unlist(Filter(Negate(is.null), corrupt_files))
#length(corrupt_files)
#file.remove(corrupt_files) # Remove corrupt files

#-#-# Add area size, clip and calculate loss of climatically suitable space #-#-#
n <- ceiling(0.75*parallel::detectCores())
sfInit(parallel=TRUE, cpus=n)
sfLibrary(raster); sfLibrary(sp); sfLibrary(dplyr); sfLibrary(tidyr)
sfLibrary(rgeos); sfLibrary(readr)
sfExport(list=c("area", "SpeciesData", "predpath", "resultspath", "model_type")) 
sfLapply(sample(names_mis), function(x){
  print(x)
  if(!file.exists(paste0(resultspath,  x, "_", model_type, "_dispersal.csv.xz"))){
    ## Future dist
    spPred <- readr::read_csv(paste0(predpath,x,"_", model_type, "_predict.csv.xz"))
    
    ## Current dist
    spOrig <- raster::raster(paste0(SpeciesData,x,"_0.5.tif"))
    coord <- round(coordinates(spOrig),4)
    values <- getValues(spOrig)
    origDist <- (as.data.frame(cbind(coord,values)))
    origDist <- subset(origDist,values==1)
    colnames(origDist) <- c("x","y","presence")
    
    # Calculate dispersal distance
    poly <- raster::rasterToPolygons(spOrig,fun=function(x){x == 1},dissolve=TRUE)
    poly <- sp::disaggregate(poly)
    largepoly <- poly[which.max(sapply(poly@polygons, function(x) x@Polygons[[1]]@area)),]
    xdist <- xmax(largepoly) - xmin(largepoly)
    ydist <- ymax(largepoly) - ymin(largepoly)
    
    # Various distances
    dist1 <- sqrt((xdist^2+ydist^2))/4
    disp1 <- rgeos::gBuffer(poly, byid=TRUE, width=dist1)
    disp1 <- crop(disp1, extent(spOrig))
    disp1 <- raster::rasterize(disp1, spOrig, field=1)
    dispDist1 <- as.data.frame(raster::rasterToPoints(disp1))
    colnames(dispDist1) <- c("x","y","dispersal1")
    
    dist2 <- sqrt((xdist^2+ydist^2))/2
    disp2 <- rgeos::gBuffer(poly, byid=TRUE, width=dist2)
    disp2 <- crop(disp2, extent(spOrig))
    disp2 <- raster::rasterize(disp2, spOrig, field=1)
    dispDist2 <- as.data.frame(rasterToPoints(disp2))
    colnames(dispDist2) <- c("x","y","dispersal2")
    
    dist3 <- sqrt((xdist^2+ydist^2))
    disp3 <- rgeos::gBuffer(poly, byid=TRUE, width=dist3)
    disp3 <- crop(disp3, extent(spOrig))
    disp3 <- raster::rasterize(disp3, spOrig, field=1)
    dispDist3 <- as.data.frame(raster::rasterToPoints(disp3))
    colnames(dispDist3) <- c("x","y","dispersal3")
    
    dist4 <- sqrt((xdist^2+ydist^2))*2
    disp4 <- rgeos::gBuffer(poly, byid=TRUE, width=dist4)
    disp4 <- crop(disp4, extent(spOrig))
    disp4 <- raster::rasterize(disp4, spOrig, field=1)
    dispDist4 <- as.data.frame(raster::rasterToPoints(disp4))
    colnames(dispDist4) <- c("x","y","dispersal4")
    
    dataList <- list(dispDist4,dispDist3,dispDist2,dispDist1)
    
    Dispersal <- Reduce(function(...) merge(..., all=T), dataList)
    
    ## Merge dist files
    CP <- dplyr::full_join(origDist, Dispersal, by=c("x","y"))
    CP <- dplyr::left_join(spPred, CP, by=c("x", "y"))
    CP$fulldisp <- 1
    
    ## Add cell size 
    CPA <- dplyr::left_join(CP, area, by=c("x","y"))
    
    ## Set 0 values to NA
    CPA[CPA == 0] <- NA
    
    # Remove rows where all values are NA
    # CPA <- tidyr::drop_na(CPA, -c("x","y", "presence", 
    # "dispersal1", "dispersal2", "dispersal3", "dispersal4", "fulldisp"))
    CPA <- CPA[apply(CPA, 1, function(y) !all(is.na(y))),]
    
    # Change column names to correct names
    colnames(CPA) <- sub("-", ".", colnames(CPA))
    
    # Remove dispersal scenario for 1995
    #CPA$EWEMBI_1995[is.na(CPA$presence)] <- 0
    
    # Save to file
    readr::write_csv(CPA, paste0(resultspath,  x, "_", model_type, 
                                 "_dispersal.csv.xz"))
    return(x)
  }
})
sfStop()
