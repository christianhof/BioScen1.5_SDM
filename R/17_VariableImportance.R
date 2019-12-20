# Variable Importance GBM

# Set file directory
#filedir <- "/bigdata/mbiber" # shinichi
#filedir <- "/scratch/mbiber/data" # darkstar
#filedir <- "/bigdata_local/mbiber" # ceremony

# Set taxa
taxa <- c("Amphibian", "Ter_Mammal", "Ter_Bird")
i <- 2

# Model type
model_type <- "GBM"

# Ecoblocking (1) or 30/70 Split (2) or Eco-Missing (3)
j <- 3

# List all files
modelpath <- c(paste0(filedir, "/", taxa[i], "_", model_type, "_plots"), 
               paste0(filedir, "/", taxa[i], "_", model_type, "_plots_30_70"),
               paste0(filedir, "/", taxa[i], "_", model_type, "_plots_30_70_MissingEco"))[j]
Modelfiles <- list.files(path=modelpath, pattern="_relative_influence.*.txt", full.names=TRUE)
head(Modelfiles)
length(Modelfiles)

# AUC Summary paths
var_inf_all <- c(paste0(filedir, "/VariableInfluence_All_", model_type, "_", taxa[i], ".csv"), 
             paste0(filedir, "/VariableInfluence_All_", model_type, "_30_70_", taxa[i], ".csv"),
             paste0(filedir, "/VariableInfluence_All_", model_type, "_30_70_MissingEco_", taxa[i], ".csv"))[j]
var_inf_sp <- c(paste0(filedir,"/VariableInfluence_Sp_", model_type, "_", taxa[i], ".csv"),
            paste0(filedir,"/VariableInfluence_Sp_", model_type, "_30_70_", taxa[i], ".csv"),
            paste0(filedir,"/VariableInfluence_Sp_", model_type, "_30_70_Missing_Eco_", taxa[i], ".csv"))[j]

# Try on one output
Model1File <- read.delim(Modelfiles[1])
head(Model1File)
str(Model1File)

#' Loop through all model output files to extract AUC and save summarized output
library(snowfall)
sfInit(parallel=TRUE, cpus=ceiling(0.9*parallel::detectCores()))
all.var.inf <- sfLapply(Modelfiles,function(modList){
  print(modList)
  
  ##Load species model
  data <- read.delim(modList)
  
  ##Extract Variables
  data <- tidyr::separate(data, var.rel.inf, into=c("row.names", "var", "rel_inf"), sep=" ")
  
  ##Get Block number 
  data$Block <- strsplit(basename(modList),split="_",fixed=TRUE)[[1]][4]
  
  ##Get PA name 
  data$PA <- strsplit(basename(modList),split="_",fixed=TRUE)[[1]][3]
  
  # Get Species name
  data$species <- paste0(strsplit(basename(modList),split="_",fixed=TRUE)[[1]][1:2],collapse="_")
  return(data)
}); sfStop()
all.var.inf <- Filter(Negate(is.null), all.var.inf)
all.var.inf <- do.call(rbind, all.var.inf)
# Change format of dataframe
all.var.inf <- tidyr::spread(all.var.inf, Block, rel_inf)
colnames(all.var.inf)<-c("row.names", "var", "PA", "Species", 
                         "Inf.1","Inf.2","Inf.3","Inf.4","Inf.5",
                         "Inf.6","Inf.7","Inf.8","Inf.9","Inf.10")
head(all.var.inf)
nrow(all.var.inf)
## Save AUC file as csv
readr::write_csv(all.var.inf, 
                 path=var_inf_all); rm(all.var.inf)

# Aggregate the different AUC values from the 10 iterations per species
all.var.inf <- read.csv(var_inf_all)
head(all.var.inf)

Infdata <- all.var.inf[,c(5:ncol(all.var.inf))]
Infdata <- rowMeans(Infdata, na.rm=TRUE)
all.var.inf$MeanInf <- Infdata
All.Sub <- all.var.inf[c(2,3,4,ncol(all.var.inf))]
All.Sub <- tidyr::spread(All.Sub, var, MeanInf)
All.Sub <- All.Sub[c(2:ncol(All.Sub))]
head(All.Sub)
nrow(All.Sub)
All.Inf.MeanPerSp <- aggregate(.~Species, data=All.Sub, mean, na.action= na.omit)
head(All.Inf.MeanPerSp)
nrow(All.Inf.MeanPerSp)

readr::write_csv(All.Inf.MeanPerSp, 
                 path=var_inf_sp)

# Merge AUC files from the different model_types (eco, 30_70, 30_70_Eco_missing)

# Set file directory
filedir <- "E:/ProcessedData" # Desktop

# Set taxa
taxa <- c("Amphibian", "Ter_Mammal", "Ter_Bird")
i <- 2

# Model type
model_type <- "GBM"

# AUC Summary paths
var_inf_all <- c(paste0(filedir, "/VariableInfluence/VariableInfluence_All_", model_type, "_", taxa[i], ".csv"), 
                 paste0(filedir, "/VariableInfluence/VariableInfluence_All_", model_type, "_30_70_", taxa[i], ".csv"),
                 paste0(filedir, "/VariableInfluence/VariableInfluence_All_", model_type, "_30_70_MissingEco_", taxa[i], ".csv"))
var_inf_sp <- c(paste0(filedir,"/VariableInfluence/VariableInfluence_Sp_", model_type, "_", taxa[i], ".csv"),
                paste0(filedir,"/VariableInfluence/VariableInfluence_Sp_", model_type, "_30_70_", taxa[i], ".csv"),
                paste0(filedir,"/VariableInfluence/VariableInfluence_Sp_", model_type, "_30_70_Missing_Eco_", taxa[i], ".csv"))

# Merge files
library(dplyr)
inf_all <- lapply(1:3, function(x){
  data <- read.csv(var_inf_all[x])
  data$type <- c("Eco", "30_70", "Missing_Eco")[x]
  return(data)
})
org_name <- setdiff(inf_all[[1]][,c(1,2)], inf_all[[2]][,c(1,2)])
org_name2 <- setdiff(inf_all[[1]][,c(1,2)], inf_all[[3]][,c(1,2)])
org_name <- intersect(org_name, org_name2)
inf_eco <- inf_all[[1]][inf_all[[1]][,"Species"] %in% org_name$Species & 
                          inf_all[[1]][,"PA"] %in% org_name$Iteration,]
org_name <- setdiff(inf_all[[2]][,c(1,2)], inf_all[[3]][,c(1,2)])
inf_30_70 <- inf_all[[2]][inf_all[[2]][,"Species"] %in% org_name$Species & 
                            inf_all[[2]][,"PA"] %in% org_name$Iteration,]
inf_all <- do.call("rbind", list(inf_eco, inf_30_70, inf_all[[3]]))
which(duplicated(inf_all[,c("var", "Species", "PA")]))

inf_sp <- lapply(1:3, function(x){
  data <- read.csv(var_inf_sp[x])
  data$type <- c("Eco", "30_70", "Missing_Eco")[x]
  return(data)
})
org_name <- setdiff(inf_sp[[1]][,1], inf_sp[[2]][,1])
org_name2 <- setdiff(inf_sp[[1]][,1], inf_sp[[3]][,1])
org_name <- intersect(org_name, org_name2)
inf_eco <- inf_sp[[1]][inf_sp[[1]][,"Species"] %in% org_name,]
org_name <- setdiff(inf_sp[[2]][,1], inf_sp[[3]][,1])
inf_30_70 <- inf_sp[[2]][inf_sp[[2]][,"Species"] %in% org_name,]
inf_sp <- do.call("rbind", list(inf_eco, inf_30_70, inf_sp[[3]]))
which(duplicated(inf_sp[,"Species"]))

# Save files and delete old files
write.csv(inf_all, paste0(filedir, "/VariableInfluence/Var_Inf_all_", taxa[i], "_", model_type, ".csv"), row.names=FALSE)
write.csv(inf_sp, paste0(filedir, "/VariableInfluence/Var_Inf_sp_", taxa[i], "_", model_type, ".csv"), row.names=FALSE)
# Delete individual files
#file.remove(auc_all_files)
#file.remove(auc_sp_files)

# Create plot of variable influence

inf_all <- lapply(taxa, function(x){
  data <- read.csv(paste0(filedir, "/VariableInfluence/Var_Inf_all_", x, "_GBM.csv"))
  data$taxa <- x
  return(data)
})
inf_all <- do.call("rbind", inf_all)
inf_all <- tidyr::gather(inf_all, Block, value, -c(row.names, var, PA, Species, type, taxa))
inf_all$var <- factor(inf_all$var, levels=c("bio4", "bio5", "bio12", "bio15", "bio18", "bio19"))
inf_all$taxa <- factor(inf_all$taxa, labels=c("Amphibians", "Birds", "Mammals"))
library(ggplot2)
ggplot(data=inf_all, aes(x=var, y=value, fill=taxa)) + geom_boxplot() +
  theme_classic() + labs(x="Variable", y="Relative Influence") + 
  scale_fill_manual(name="Taxa", values=c(rgb(0.3,0.5,0.1,1), rgb(0.8,0.1,0.1,1), 
                                        rgb(0.2,0.2,0.6,1)))
ggsave("figures/Rel_Influence_GBM_All.png", width=8, height=5, dpi=100)

inf_sp <- lapply(taxa, function(x){
  data <- read.csv(paste0(filedir, "/VariableInfluence/Var_Inf_sp_", x, "_GBM.csv"))
  data$taxa <- x
  return(data)
})
inf_sp <- do.call(plyr::rbind.fill, inf_sp)
inf_sp <- tidyr::gather(inf_sp, var, value, -c(Species, type, taxa))
inf_sp$var <- factor(inf_sp$var, levels=c("bio4", "bio5", "bio12", "bio15", "bio18", "bio19"))
inf_sp$taxa <- factor(inf_sp$taxa, labels=c("Amphibians", "Birds", "Mammals"))
ggplot(data=inf_sp, aes(x=var, y=value, fill=taxa)) + geom_boxplot() +
  theme_classic() + labs(x="Variable", y="Relative Influence") + 
  scale_fill_manual(name="Taxa", values=c(rgb(0.3,0.5,0.1,1), rgb(0.8,0.1,0.1,1), 
                                          rgb(0.2,0.2,0.6,1)))
ggsave("figures/Rel_Influence_GBM_Sp.png", width=8, height=5, dpi=100)
