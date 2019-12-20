#-#-# Uncertainty around projected distributions #-#-#

library(ggplot2)
library(dplyr)
library(tidyr)

# Set file directory
filedir <- "/bigdata_local/mbiber/" # ceremony
#filedir <- "J:/"

# Set working directory
setwd("/home/mbiber/BioScen1.5_1") #ceremony

## GAM file path
filepaths <- c("Ter_Mammal_GAM_results_climate/", "Ter_Bird_GAM_results_climate/", 
               "Amphibian_GAM_results_climate/", "Ter_Mammal_GBM_results_climate/", 
               "Ter_Bird_GBM_results_climate/", "Amphibian_GBM_results_climate/")

## Líst all files 
allFiles <- sapply(filepaths,function(x){list.files(paste0(filedir, x), full.names=T)})
allFiles <- do.call(c, allFiles)
length(allFiles)

## Extract plot data
library(snowfall)
sfInit(parallel=TRUE, cpus=ceiling(0.75*parallel::detectCores()))
sfLibrary(dplyr); sfLibrary(tidyr)
sfExport(list=c("filedir")) 

plotData <- sfLapply(allFiles, function(x){
  print(x)
  spdata <- read.csv(x)
  if(nrow(spdata) > 0){
  head(spdata)
  #spdata <- na.omit(spdata)
  
  spdata$species <- paste0(strsplit(basename(x), split="_")[[1]][1:2], collapse="_")
  spdata$model_type <- strsplit(basename(x), split="_")[[1]][3]
  spdata$taxa <- strsplit(sub("Ter_", "", sub(filedir, "", dirname(x))),split="_")[[1]][1]
  
  # Subset spdata
  spdata_all <- spdata %>% select(x,y,presence,dispersal, species, model_type, taxa, 
                              EWEMBI_1995:MIROC5_rcp60_2080)
  
  #Calculate mean change among RCP, Year, Model algorithm and GCM
  spdata <- spdata_all %>% group_by(species, model_type, taxa) %>% 
    summarise_at(vars(`EWEMBI_1995`:`MIROC5_rcp60_2080`), mean, na.rm=T) %>%
    mutate_at(vars(`GFDL.ESM2M_rcp26_2050`:`MIROC5_rcp60_2080`), funs(. - EWEMBI_1995)) %>% 
    select(-EWEMBI_1995)
  spdata$disp <- 1
  if(nrow(spdata_all %>% filter(presence == 1)) > 0){
    spdata_nodisp <- spdata_all %>% filter(presence == 1) %>% 
      group_by(species, model_type, taxa) %>% 
      summarise_at(vars(`EWEMBI_1995`:`MIROC5_rcp60_2080`), mean, na.rm=T) %>%
      mutate_at(vars(`GFDL.ESM2M_rcp26_2050`:`MIROC5_rcp60_2080`), funs(. - EWEMBI_1995)) %>% 
      select(-EWEMBI_1995)
    spdata_nodisp$disp <- 0
    spdata <- rbind(spdata, spdata_nodisp)    
  }
  return(spdata)
  }
})
meanprobabilitychange <- data.table::rbindlist(plotData)
head(meanprobabilitychange)
str(meanprobabilitychange)
write.csv2(meanprobabilitychange, paste0(filedir, "/meanprobabilitychange.csv"), row.names=F)

#-#-# Plot the mean change in probabilities
library(dplyr); library(tidyr); library(ggplot2)
meanprobabilitychange <- read.csv2("data/meanprobabilitychange.csv")
#meanprobabilitychange$taxa <- sub("/", "", meanprobabilitychange$taxa)
meanprobability <- tidyr::gather(meanprobabilitychange, var, value, -c(species, model_type, taxa, disp))
meanprobability <- tidyr::separate(meanprobability, var, into=c("GCM", "rcp", "time"), sep="_")
meanprobability <- tidyr::unite(meanprobability, col="ind", GCM, model_type)

lapply(c("disp", "nodisp"), function(disp){
  lapply(c(2050, 2080), function(year){
    if(disp == "nodisp"){meanchange_sub <- meanprobability %>% filter(disp == 0)} else{
      meanchange_sub <- meanprobability %>% filter(disp == 1)
    }
    meanchange_sub <- meanchange_sub %>% filter(time == year)
    meanchange_sub$rcp <- factor(meanchange_sub$rcp, labels=c("RCP2.6", "RCP6.0"))
    meanchange_sub$taxa <- factor(meanchange_sub$taxa, labels=c("Amphibians", "Birds", "Mammals"))
    meanchange_sub$ind <- factor(meanchange_sub$ind, labels=c("GFDL-ESM2M/GAM", "GFDL-ESM2M/GBM", 
                                                              "HadGEM2-ES/GAM", "HadGEM2-ES/GBM",
                                                              "IPSL-CM5A-LR/GAM", "IPSL-CM5A-LR/GBM",
                                                              "MIROC5/GAM", "MIROC5/GBM"))
    ggplot(meanchange_sub, aes(x=value)) + 
      geom_freqpoly(aes(colour=ind), binwidth=0.05) + 
      facet_grid(taxa ~ rcp, scales="free") + 
      scale_colour_manual(name= "Different predictions", 
                          values = c("red","mediumvioletred","orange","yellow","dodgerblue4",
                                     "cornflowerblue","springgreen4","limegreen")) +
      geom_vline(xintercept = 0) +
      scale_x_continuous(name="Mean change in probability", limits=c(-1,0.85), 
                         expand=c(0.01,0.01)) +
      scale_y_continuous(name="Number of species", limits=c(0,NA), expand=c(0.01,0)) + 
      theme_classic()+
      theme(axis.title.y = element_text(size = 12), 
            axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(color="black", size=8), 
            axis.text.y = element_text(color="black", size=8),
            strip.placement = "outside",
            strip.text.x = element_text(size=10, face="bold"),
            strip.text.y = element_text(size=12, face="bold"),
            strip.background= element_blank()) + 
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)
    ggsave(paste0("figures/uncertainty_probchange_", year, "_", disp, ".png"), 
           width=8, height=6, dpi=600, bg="transparent")
  })
})
