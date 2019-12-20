##################################################

#-#-# Climate data stability #-#-#
## - using PCA of climate variables and Euclidean distances between PC of the curent and future climate

rm(list=ls(all=TRUE))

## Load libraries
library(psych)
library(lattice)
library(dplyr)
library(ggplot2)
library(tidyr)

## Get dataframes current and future
## Current climate 
bioclim_1995 <- get(load(list.files(path="/home/mbiber/Documents/GitHub/rISIMIP/data/", 
                                    pattern="bioclim_EWEMBI.*\\.rda", 
                                    full.names=T, recursive=T)))[,c("x", "y", "bio4", "bio5",  
                                                                   "bio12", "bio15", "bio18", 
                                                                   "bio19")]

## PC current climate data
pc_bird_mam <- prcomp(bioclim_1995[,c("bio4", "bio5", "bio12", "bio15")], center= TRUE, scale = TRUE)
pc_amphi <- prcomp(bioclim_1995[,c("bio4", "bio5", "bio18", "bio19")], center= TRUE, scale = TRUE)

#plot(pc_bird_mam) #Which principal component captures most of the variance
#plot(pc_bird_mam,type="l") #Ellbow method to decide on number of components - Bend at 3 
#biplot(pc_bird_mam) #Shows bearing of the input - little correlation between components
#dim(bioclim_1995[3:6]) # 67420 observations of 4 variables
#pairs.panels(bioclim_1995[3:6], gap = 0)

## Extract the current PC data
#str(pc)

climdata_bird <- cbind(bioclim_1995,pc_bird_mam$x[,c("PC1", "PC2")])
climdata_bird$taxa <- "Ter_Birds"
climdata_mam <- cbind(bioclim_1995,pc_bird_mam$x[,c("PC1", "PC2")])
climdata_mam$taxa <- "Ter_Mammals"
climdata_amphi <- cbind(bioclim_1995, pc_amphi$x[,c("PC1", "PC2")])
climdata_amphi$taxa <- "Amphibians"
climdata <- rbind(climdata_bird, climdata_mam, climdata_amphi)
  
## Plot current PC data
ggplot(climdata, aes(x=PC1,y=PC2))+ facet_wrap(~ taxa) + 
  stat_ellipse(geom="polygon", col = "black", alpha=0.5)+
  geom_point(shape=1)

## Return final data
CurrentData <- climdata[c("x", "y","taxa", "PC1","PC2")] 
CurrentData$year <- 1995
head(CurrentData)

## Future climate
bioclim_fut <- c(list.files(path="/home/mbiber/Documents/GitHub/rISIMIP/data/", 
                            pattern="bioclim_.*2050.*\\.rda", full.names=T, recursive=T), 
                 list.files(path="/home/mbiber/Documents/GitHub/rISIMIP/data/", 
                            pattern="bioclim_.*2080.*\\.rda", full.names=T, recursive=T))
bioclim_fut <- lapply(bioclim_fut, function(x){
  data <- get(load(x))
  data$year <- strsplit(strsplit(basename(x), split="_")[[1]][4], split="[.]")[[1]][1]
  data$model <- strsplit(basename(x), split="_")[[1]][2]
  data$scenario <- strsplit(basename(x), split="_")[[1]][3]
  return(data)
})
bioclim_fut <- do.call("rbind", bioclim_fut)

# Calculate mean across GCMs
#bioclim_fut <- tidyr::gather(bioclim_fut, var, value, -c(x,y,model, scenario, year))
#bioclim_fut <- bioclim_fut %>% group_by(x,y,scenario,year,var) %>% 
#  summarise(mean = mean(value, na.rm=TRUE))
#bioclim_fut <- tidyr::spread(bioclim_fut, var, mean)
#head(bioclim_fut)

## Run PCA separated by scenario and year
FutureData <- lapply(c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5"), function(model){
 yeardata <- lapply(c(2050,2080), function(time){
  climdata <- lapply(c("rcp26", "rcp60"), function(rcp){
    bio_sub_bird_mam <- bioclim_fut[bioclim_fut$year== time & bioclim_fut$model== model & 
                                           bioclim_fut$scenario == rcp, 
                                         c("x", "y", "bio4","bio5","bio12","bio15")]
    pc_bird_mam <- prcomp(bio_sub_bird_mam[, c("bio4","bio5","bio12","bio15")], 
                          center= TRUE, scale = TRUE)
    bio_sub_amphi <- bioclim_fut[bioclim_fut$year== time & bioclim_fut$model== model & 
                                         bioclim_fut$scenario == rcp,
                                       c("x", "y", "bio4","bio5", "bio18","bio19")]
    pc_amphi <- prcomp(bio_sub_amphi[,c("bio4","bio5","bio18","bio19")], 
                       center= TRUE, scale = TRUE)
    climdata_bird <- cbind(bio_sub_bird_mam,pc_bird_mam$x[,c("PC1", "PC2")])
    climdata_bird$taxa <- "Ter_Birds"
    climdata_mam <- cbind(bio_sub_bird_mam,pc_bird_mam$x[,c("PC1", "PC2")])
    climdata_mam$taxa <- "Ter_Mammals"
    climdata_amphi <- cbind(bio_sub_amphi, pc_amphi$x[,c("PC1", "PC2")])
    climdata_amphi$taxa <- "Amphibians"
    climdata <- plyr::rbind.fill(climdata_bird, climdata_mam, climdata_amphi)
    climdata$scenario <- rcp
    return(climdata)
  })
  climdata <- do.call("rbind",climdata)
  climdata$year <- time
  return(climdata)
})
  yeardata <- do.call("rbind", yeardata)
  yeardata$model <- model
  return(yeardata)
})
FutureData <- do.call("rbind", FutureData)
FutureData <- FutureData[c("x", "y","taxa", "PC1","PC2", "year", "scenario", "model")] 

## Plot future PC data
#ggplot(FutureData, aes(x=PC1,y=PC2)) + facet_grid(year + scenario + model ~ taxa) + 
#  stat_ellipse(geom="polygon", col = "black", alpha=0.5)+
#  geom_point(shape=1)

## Merge current and future
CurrentData1 <- CurrentData
CurrentData1$scenario <- "rcp26"
CurrentData2 <- CurrentData
CurrentData2$scenario <- "rcp60"
CurrentData <- rbind(CurrentData1, CurrentData2)
CurrentData <- lapply(c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5"), function(model){
  data <- CurrentData
  data$model <- model
  return(data)
})
CurrentData <- bind_rows(CurrentData)
bioclim_all <- rbind(CurrentData, FutureData)
head(bioclim_all)

## Subset for 2050 and 2080
bioclim_2050 <- bioclim_all %>% filter(year != 2080)
bioclim_2080 <- bioclim_all %>% filter(year != 2050)

edist_2050 <- bioclim_2050 %>% ungroup() %>% select(-year) %>% 
  group_by(x,y, scenario, taxa, model) %>% 
  dplyr::summarise(dist=dist(cbind(PC1, PC2), method="euclidean"))
edist_2050$year <- 2050
edist_2080 <- bioclim_2080 %>% group_by(x,y, scenario, taxa, model) %>% 
  summarise(dist=dist(cbind(PC1, PC2), method="euclidean"))
edist_2080$year <- 2080
edist <- rbind(edist_2050, edist_2080)

# Save to file
readr::write_csv(edist, "data/climate_stability.csv.xz")

# Plot climate stability
library(ggplot2)
ggplot() + geom_tile(data=edist, aes(x=x, y=y, fill=dist)) + 
  facet_grid(year + model ~taxa) + scale_fill_gradientn(colours=rainbow(255))

##################################################

# Calculate climate top25

edist <- readr::read_csv("data/climate_stability.csv.xz")

#define climate change thresholds
summary(edist$dist)
hist(edist$dist)
(thres <- edist %>% group_by(year, taxa) %>% summarise(thres = quantile(dist, probs=0.75)))

# Merge threshold with data
edist <- left_join(edist, thres)

# Subset data by threshold
edist <- edist %>% filter(dist > thres)
hist(edist$dist)

# Plot
ggplot() + geom_tile(data=edist, aes(x=x, y=y, fill=dist)) + 
  facet_grid(year~taxa) + scale_fill_gradientn(colours=rainbow(255))

# Save subset to file
readr::write_csv(edist, "data/top_climate_stability.csv.xz")

##################################################