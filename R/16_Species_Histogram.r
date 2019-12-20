########################################

## Create climate and landuse summaries

library(snowfall)
library(dplyr)
library(tidyr)

# Set file directory
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/bigdata_local/mbiber" # ceremony - Mammals
filedir <- "/home/mbiber/data" # ceremony - Birds

#filedir <- "E:/ProcessedData"

# Set working directory
#setwd("/scratch/home/mbiber/GitHub/BioScen1.5_1") #shinichi
#setwd("/home/mbiber/BioScen1.5_1") #ceremony

# Set taxa
k <- 3
i <- c("Amphibian", "Mammal", "Bird")[k]

# Get available files
files_GAM <- list.files(paste0(filedir, "/", i, "_GAM_results_climate/"), pattern=".csv",
                        full.names=T)
length(files_GAM)
files_GBM <- list.files(paste0(filedir, "/", i, "_GBM_results_climate/"), pattern=".csv",
                        full.names=T)
length(files_GBM)

# Only take species that have been modelled
names_gam <- readRDS("data/names_gam.rds")[[k]]
names_gbm <- readRDS("data/names_gbm.rds")[[k]]
names_mod <- names_gam[names_gbm %in% names_gam]

names <- lapply(files_GAM , function(x) paste0(strsplit(basename(x), "_")[[1]][1:2], 
                                               collapse="_"))
files_GAM <- files_GAM[names %in% names_mod]
length(files_GAM)
names <- lapply(files_GBM, function(x) paste0(strsplit(basename(x), "_")[[1]][1:2], 
                                              collapse="_"))
files_GBM <- files_GBM[names %in% names_mod]
length(files_GBM)

# Remove files with low AUC

# Get low AUC species
lowAUC <- read.csv2("data/species_lowAUC.csv")

# Subset by taxa
lowAUC <- lowAUC[lowAUC$taxa == i,]
nrow(lowAUC)

# Extract names and subset files by names that have low AUC
names <- lapply(files_GAM , function(x) paste0(strsplit(basename(x), "_")[[1]][1:2], 
                                               collapse="_"))
files_GAM <- files_GAM[!names %in% lowAUC$Species[lowAUC$model_type == "GAM"]]
length(files_GAM)
names <- lapply(files_GBM, function(x) paste0(strsplit(basename(x), "_ma_")[[1]][1:2], 
                                              collapse="_"))
files_GBM <- files_GBM[!names %in% lowAUC$Species[lowAUC$model_type == "GBM"]]
length(files_GBM)

# Extract all unique species
spname_GAM <- lapply(files_GAM, function(x){
  paste(strsplit(basename(x),split="_")[[1]][1:2],collapse="_")
})
spname_GBM <- lapply(files_GBM, function(x){
  paste(strsplit(basename(x),split="_")[[1]][1:2],collapse="_")
})
spname <- unique(c(spname_GAM, spname_GBM))
length(spname)

# Drop two bird species
spname <- spname[!spname %in% "Arborophila_tonkinensis"]
spname <- spname[!spname %in% "Phylloscartes_parkeri"]
spname <- spname[!spname %in% "Pipra_fasciicauda"]
spname <- spname[!spname %in% "Ploceus_nicolli"]
spname <- spname[!spname %in% "Xiphorhynchus_ocellatus"]

# Read land_use data
landuse_all <- read.csv("data/landuse_all.csv.xz")
# Do not change to readr::read_csv, does not work somehow.

# Change land-use categories
landuse_all$biofuel_cropland <- landuse_all$biofuel_cropland_irrigated + landuse_all$biofuel_cropland_rainfed
landuse_all$biofuel_cropland[is.na(landuse_all$biofuel_cropland)] <- 0
landuse_all$cropland <- landuse_all$cropland_rainfed + landuse_all$cropland_irrigated

# Change format of land-use data
landuse_all %<>% select(-c(cropland_irrigated, cropland_rainfed,
                           biofuel_cropland_irrigated, 
                           biofuel_cropland_rainfed)) %>% 
  tidyr::gather(var, value, -c(x,y,year, scenario, model)) %>% tidyr::spread(year, value)

# Calculate delta landuse
delta_landuse <- landuse_all %>% 
  mutate_at(vars(`2050`:`2080`), funs(. - `1995`)) %>% 
  dplyr::select(c("x", "y", "scenario", "var", "model", "2050", "2080")) %>% 
  tidyr::gather(year, value, -c(x, y, scenario, var, model))

# Calculate ensemble mean
delta_mean <- delta_landuse %>% group_by(x,y,scenario,var,year) %>% 
  dplyr::summarise(mean=mean(value, na.rm=TRUE))

# Define overall landuse threshhold
(threshold <- delta_landuse %>% group_by(var, year) %>% filter(value > 0) %>% 
    summarise(threshold = quantile(value, probs=0.75)))

# Define landuse threshold for each GCM and algorithm
#(threshold <- delta_landuse %>% group_by(var, year, model) %>% filter(value > 0) %>% 
#    summarise(threshold = quantile(value, probs=0.75)))

# Define mean land use threshold
#(threshold <- delta_mean %>% group_by(var, year) %>% filter(mean > 0) %>% 
#    summarise(threshold = quantile(mean, probs=0.75)))

# Subset data by threshold
landuse_change <- left_join(delta_landuse, threshold, by=c("var", "year")) %>% 
  group_by(x,y,scenario, model, var, year) %>% 
  filter(value > threshold) %>% select(-threshold) %>% mutate(value = 1)

# Adapt format of landuse data (from long to wide)
landuse_change <- tidyr::spread(landuse_change, var, value)
landuse_change$year <- as.character(landuse_change$year)

# Turn numbers into true or false
landuse_change$biofuel_cropland <- landuse_change$biofuel_cropland > 0
landuse_change$cropland <- landuse_change$cropland > 0
landuse_change$pastures <- landuse_change$pastures > 0

landuse_change$GCM <- landuse_change$model
landuse_change$GCM <- factor(landuse_change$GCM, labels=c("GFDL.ESM2M", "HadGEM2.ES", "IPSL.CM5A.LR", "MIROC5"))
landuse_change2 <- landuse_change
landuse_change$model <- "GAM"
landuse_change2$model <- "GBM"
landuse_change <- bind_rows(landuse_change, landuse_change2)

sfInit(parallel=TRUE, cpus=ceiling(0.85*parallel::detectCores()))
sfLibrary(dplyr); sfLibrary(tidyr)
sfExport(list=c("files_GAM", "files_GBM", "landuse_change")) 

source("R/clim_lu_impact_func.R")
ChangeSums <- sfLapply(spname, clim_lu_impact)
ChangeSums <- Filter(Negate(is.null), ChangeSums)
AlldataLUCS <- data.table::rbindlist(ChangeSums, fill=T)
head(AlldataLUCS)
readr::write_csv(AlldataLUCS, paste0(filedir, "/", i, "_ClimateLU_Impact.csv.xz"))
sfStop()

########################################

## Create frequency plot (or histogram) and bar graph

# Load packages
library(dplyr); library(tidyr); library(ggplot2); library(magrittr)

# Set year, dispersal and threshold
time <- 2080
disp <- "nodisp"
threshold <- 0 

# Read data
#a <- read.csv(paste0(filedir, "/Amphibian_ClimateLU_Impact.csv"), as.is=T)
a <- read.csv(paste0("data/Amphibian_ClimateLU_Impact.csv.xz"), as.is=T)
a$taxa <- "Amphibians"
b <- read.csv(paste0("data/Bird_ClimateLU_Impact.csv.xz"), as.is=T)
b$taxa <- "Birds"
m <- read.csv(paste0("data/Mammal_ClimateLU_Impact.csv.xz"), as.is=T)
m$taxa <- "Mammals"
d <- dplyr::bind_rows(list(a,b,m)); rm(a,b,m)

# Subset by threshold and year
d %<>% filter(thres == threshold) %>% filter(year == time)

# Subset by dispersal
if(disp == "disp"){
  d %<>% filter(disp == 1)
} else{
  d %<>% filter(disp == 0)
}

# Identify number of species considered
length(unique(d$species))

# Calculate proportion of area
d %<>% ungroup %>% group_by(scenario, year, model, GCM, species, taxa) %>% 
  dplyr::mutate(total=sum(sum)) %>% mutate(prop=sum/total)

d$threat <- factor(d$threat, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                   labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
                            "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
                            "CC+BC+PA", "CC+CR+PA", "BC+CR+PA", "CC+BC+CR+PA"))

mean_threat <- d %>% group_by(scenario, year, taxa, species, threat) %>% 
  summarise(tp=mean(prop)) %>% tidyr::spread(threat, tp)
readr::write_csv(mean_threat, "mean_threat_nodisp.csv")

mean_threat <- read.csv("mean_threat_nodisp.csv")

mean_threat_species <- mean_threat %>% 
  filter(species %in% c("Perdix_perdix", "Alauda_arvensis",
                        "Vanellus_vanellus", "Emberiza_citrinella"))
readr::write_csv(mean_threat_species, "mean_threat_nodisp_species.csv")

# Combine data with threat possibilities
df <- expand.grid(taxa=c("Amphibians", "Birds", "Mammals"),
                  scenario = c("rcp26","rcp60"), year=time,
                  model = c("GAM", "GBM"), GCM=c("MIROC5", "GFDL.ESM2M", "HadGEM2.ES", "IPSL.CM5A.LR"),
                  threat=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,NA))
d <- full_join(d,df)

# Subset data by rcp and year
d$scenario <- factor(d$scenario, labels=c("RCP2.6", "RCP6.0"))
d$perc <- round(d$prop*10)/10
d$perc[d$perc == 0] <- 0.1
d$perc <- d$perc - 0.05
d$threat <- factor(d$threat, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                   labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
                            "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
                            "CC+BC+PA", "CC+CR+PA", "BC+CR+PA", "CC+BC+CR+PA"))

# Calculate number of species that are affected by a combination of impacts.
d %>% filter(threat %in% c("CC", "BC", "CC+BC")) %>% 
  group_by(scenario, model, GCM) %>% 
  summarise(count=length(unique(species))) %>% ungroup() %>% 
  group_by(scenario) %>%
  dplyr::summarise(mean=mean(count,na.rm=T), sd=sd(count,na.rm=T))

d_new <- d %>% group_by(taxa, scenario, threat, perc, model, GCM) %>% 
  dplyr::summarise(count = n()) %>% tidyr::drop_na() %>% 
  group_by(taxa, scenario, threat, perc) %>% 
  dplyr::summarise(mean=mean(count,na.rm=T), sd=sd(count,na.rm=T))

# Check number of species that are affected by Biofuel 
# by more than 10% of their range
d %>% filter(threat == "CC") %>% 
  group_by(scenario, model, GCM) %>% 
  summarise(count=sum(count)) %>% ungroup() %>% group_by(scenario) %>%
  dplyr::summarise(mean=mean(count,na.rm=T), sd=sd(count,na.rm=T))
d %>% filter(threat == "BC") %>% 
  group_by(scenario, model, GCM) %>% 
  summarise(count=sum(count)) %>% ungroup() %>% group_by(scenario) %>%
  dplyr::summarise(mean=mean(count,na.rm=T), sd=sd(count,na.rm=T))

d %>% filter(threat == "BC") %>% 
  filter(perc >= 0.15) %>% group_by(scenario, model, GCM) %>% 
  summarise(count=sum(count)) %>% ungroup() %>% group_by(scenario) %>%
  dplyr::summarise(mean=mean(count,na.rm=T), sd=sd(count,na.rm=T))

d %>% filter(threat == "CC") %>% 
  filter(perc >= 0.5) %>% group_by(scenario, model, GCM) %>% 
  summarise(count=sum(count)) %>% ungroup() %>% group_by(scenario) %>%
  dplyr::summarise(mean=mean(count,na.rm=T), sd=sd(count,na.rm=T))

# Important for right format of facetting!!!
d_new$threat <- factor(d_new$threat,levels=c('CC','BC','CR','PA','',' ',
                                             'CC+BC', 'CC+CR', 'CC+PA', 'BC+CR', 'BC+PA', 'CR+PA', 
                                             'CC+BC+CR', 'CC+BC+PA', 'CC+CR+PA', 'BC+CR+PA', '  ', 
                                             'CC+BC+CR+PA'))
# Calc error bar position
#library(plyr)
#d_new %<>% arrange(desc(taxa))
#d_new <- ddply(d_new,.(scenario,threat,perc),transform,ystart = cumsum(mean)-sd,yend = cumsum(mean) + sd)

p <- d_new %>% ggplot(aes(x=perc, y=ifelse(scenario == "RCP6.0",-1, 1)*mean, fill=taxa)) + 
  geom_col(position="stack", width=0.1) + 
  #geom_errorbar(aes(ymin=ifelse(scenario == "RCP6.0",-1, 1)*ystart, 
  #                  ymax=ifelse(scenario == "RCP6.0",-1, 1)*yend), width=0.05, colour="black") + 
  facet_wrap(~ threat, ncol=6,drop=F) + 
  #scale_fill_manual(name="Scenario", values=c("#D55E00", "#0072B2")) + 
  scale_fill_grey(name="Taxa") +
  theme_classic() + theme(strip.placement = "outside",
                          panel.spacing = unit(1, "lines"),
                          title = element_text(size=14, face="bold"),
                          plot.title = element_text(hjust = 0.5),
                          legend.text = element_text(size=16),
                          legend.title = element_text(size=16, face="bold"),
                          axis.text = element_text(size=13),
                          axis.title = element_text(size=16),
                          strip.text = element_text(size=13, face="bold"),
                          strip.background= element_blank(),
                          panel.background = element_rect(fill = "transparent", colour=NA),
                          plot.background = element_rect(fill = "transparent"),
                          legend.background = element_rect(fill = "transparent", colour=NA),
                          legend.position = c(0.85,0.9),
                          legend.box.background = element_rect(fill = "transparent", colour=NA)) + 
  scale_x_continuous(name="Proportion of distribution", labels=c(0,0.2,0.4,0.6,0.8,1),
                     breaks=seq(0,1,by=0.2), limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(name="Number of species", labels=c(9000, 6000, 3000, 0, 3000, 6000, 9000),
                     breaks=c(-9000, -6000, -3000, 0, 3000, 6000, 9000), 
                     limits=c(-9500, 9500), expand=c(0,0)) + 
  annotate("text", x=0.5, y = 5000, label = "RCP2.6", size=5) + 
  annotate("text", x=0.5, y= -5000, label= "RCP6.0", size=5) + 
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-1-3", "panel-3-3", "panel-4-3")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
## move axis closer to panel
g$layout[g$layout$name == "axis-b-5-3", c("t", "b")] = c(17, 17)

# Save plot
if(time == 2080 & disp == "nodisp" & threshold==0){
  png(paste0("figures/Figure_3.eps"),
      width=10, height=8, res=1200, units="in", bg="transparent")
  grid::grid.draw(g)
  dev.off()
  pdf(paste0("figures/Figure_3.pdf"), width=10, height=8, bg="transparent")
  grid::grid.draw(g)
  dev.off()
} else{
  png(paste0("figures/hist_clim_lu_", time, "_", disp, "_", threshold, ".png"),
      width=10, height=8, res=1200, units="in", bg="transparent")
  grid::grid.draw(g)
  dev.off()
}

