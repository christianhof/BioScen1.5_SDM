
## Get AUC stats for each model

rm(list=ls())
library(snowfall)
library(plyr)
library(data.table)

# Set file directory
#filedir <- "/scratch/home/mbiber/data/" # shinichi
#filedir <- "/bigdata_local/mbiber/" # ceremony - mammals
filedir <- "/home/mbiber/data" # ceremony - birds

#filedir <- "E:/" # Desktop

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 3

# Model type
model_type <- c("GAM", "GBM", "MaxEnt", "RF")[2]

# List all files
modelpath <- paste0(filedir, "/", taxa[i], "_", model_type, "_Output")
Modelfiles <- list.files(path=modelpath, full.names=TRUE)
length(Modelfiles)

# Try on one output
Model1File <- get(load(Modelfiles[[1]]))
length(Model1File)

## Look at AUC in first output 
Model1File[[1]][[1]]$AUC

#' Loop through all model output files to extract AUC and save summarized output
#sfInit(parallel=TRUE, cpus=ceiling(0.25*parallel::detectCores()))
#sfExport(list=c("model_type", "taxa", "i")) 

#Turn warning into error - if the model does not convert the code should stop rather than giving a warning
options(warn=2)

#510
all.sp.auc <- lapply(Modelfiles, function(modList){
  
  print(modList)
  
  # Get Species name
  spName <- paste0(strsplit(basename(modList),split="_",fixed=TRUE)[[1]][1:2],collapse="_")
  
  ##Load species model
  mtry <- try(data.list <- get(load(modList)))
  
  if(!inherits(mtry, "try-error")) {
    
    ##Extract AUC values
    AUC <- lapply(1:length(data.list), function(x){
      AUC <- sapply(1:length(data.list[[x]]), function(y){
        data.list[[x]][[y]]$AUC
      })
      AUC <- data.frame(AUC)
      AUC$PA <- paste0("PA", x)
      AUC$Block <- 1:nrow(AUC)
      return(AUC)
    })
    AUC <- Filter(Negate(is.null), AUC)
    AUC_try <- try(AUC <- do.call("rbind", AUC))
    if(!inherits(AUC_try, "try-error")) {
      AUC$Species <- spName
      AUC$taxa <- taxa[i]
      AUC$model_type <- model_type
      return(AUC)
    } else{
      return(NULL)
    }
  }
}); sfStop()
all.sp.auc <- Filter(Negate(is.null), all.sp.auc)
all.sp.auc <- data.table::rbindlist(all.sp.auc, fill=TRUE)
head(all.sp.auc)
nrow(all.sp.auc)
## Save AUC file as csv
readr::write_csv(all.sp.auc, path=
                   paste0(filedir, "/AUCvalues_All_", model_type, "_", taxa[i], ".csv")); rm(all.sp.auc)

########################################

# Create boxplot of AUC for each species and model algorithm

library(ggplot2)

# Set file directory
filedir <- "E:/ProcessedData" # Desktop

# Read data
AUC_data <- lapply(c("GAM", "GBM", "MaxEnt", "RF"), function(model_type){
  lapply(c("Amphibian", "Mammal", "Bird"), function(taxa){
    readr::read_csv(paste0(filedir, "/AUCvalues_All_", model_type, "_", taxa, ".csv"))
  })
})
AUC_data <- do.call(rbind, AUC_data)

# Aggregate the different AUC values from the 10 iterations per species
library(dplyr)
AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
  summarise(mean = mean(AUC, na.rm=T))
head(AUC_sum)

#Subset by AUC > 0.7
AUC_sum_high <- subset(AUC_sum,AUC >= 0.7)

# Make plot
ggplot(AUC_sum_high, aes(x=taxa,y=AUC,fill=model_type))+
  geom_boxplot(position=position_dodge(0.85))+
  scale_x_discrete(name = "", labels=c("Amphibians", "Birds", "Mammals")) +
  scale_fill_manual(name = "Model type", values = c("red", "blue"))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        panel.background = element_blank(),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12))
ggsave("figures/boxplot_AUC.png", width = 8, height = 6, dpi = 1200) 

########################################

# Identify species with low AUC value

# Set file directory
#filedir <- "E:/ProcessedData/AUC"

# Identify species with low AUC (smaller than 0.7)
species_lowAUC <- lapply(1:3, function(i){
  data <- lapply(1:2, function(j){
    taxa <- c("Amphibian", "Mammal", "Bird")[i]
    model_type <- c("GAM", "GBM")[j]
    data <- read.csv(paste0(filedir, "/AUC/AUC_sp_", taxa, "_", model_type, ".csv"))
    data <- dplyr::filter(data, AUC < 0.7)
    data$model_type <- model_type
    return(data)
  })
  data <- data.table::rbindlist(data)
  data$taxa <-  c("Amphibian", "Mammal", "Bird")[i]
  return(data)
})
species_lowAUC <- data.table::rbindlist(species_lowAUC)

# Save data to csv file
write.csv2(species_lowAUC, "data/species_lowAUC.csv", row.names=F)

########################################

# Remove climate files with a low AUC and move to new folder

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 3

# Model type
model_type <- c("GAM", "GBM", "RF", "MaxEnt")[4]

# Low AUC species
lowAUC <- read.csv2("data/species_lowAUC.csv")

# Subset AUC file by taxa and model_type
lowAUC <- lowAUC[lowAUC$taxa == taxa[i] | lowAUC$model_type == model_type,]

# Get all files
files <- list.files(paste0(filedir, "/", taxa, "_", model_type, "_predictions"),
                    full.names=T)

# Extract names and subset files by names that have low AUC
names <- lapply(files, function(x) paste0(strsplit(basename(x), "_")[[1]][1:2], collapse="_"))
files_sub <- files[names %in% lowAUC$Species]

# Move files
file.copy(files_sub, sub(paste0(taxa[i], "_", model_type, "_predictions"), 
                         "Low_AUC_predictions", files_sub))
file.remove(files_sub)
