library(ggplot2)
library(dplyr)
library(tidyr)

#-#-# Overlap small range species and highest change cells #-#-#

# Read climate data
climate <- readr::read_csv("data/top_climate_stability.csv.xz")
#ggplot() + geom_raster(data=climate, aes(x=x, y=y, fill=dist)) + 
#  facet_grid(scenario + year + model ~ taxa)
summary(climate[climate$taxa == "Amphibians" & climate$year == 2050,])
summary(climate[climate$taxa == "Ter_Birds" & climate$year == 2080,])

# Read landuse data
landuse_all <- read.csv("data/landuse_all.csv.xz")
# Do not use readr::read_csv here, gives error

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
rm(landuse_all)

# Define overall landuse threshhold
(threshold <- delta_landuse %>% group_by(var, year) %>% filter(value > 0) %>% 
    summarise(threshold = quantile(value, probs=0.75)))

# Subset data by threshold
landuse_change <- left_join(delta_landuse, threshold) %>% group_by(x,y,scenario, model, var, year) %>% 
  filter(value > threshold) %>% select(-threshold) %>% mutate(value = 1)
rm(delta_landuse)

# Adapt format of landuse data (from long to wide)
landuse_change <- tidyr::spread(landuse_change, var, value)
landuse_change$year <- as.character(landuse_change$year)

# Turn numbers into true or false
landuse_change$biofuel_cropland <- landuse_change$biofuel_cropland > 0
landuse_change$cropland <- landuse_change$cropland > 0
landuse_change$pastures <- landuse_change$pastures > 0

# Read landuse data
top_landuse1 <- landuse_change
top_landuse1$taxa <- "Amphibians"
top_landuse2 <- landuse_change
top_landuse2$taxa <- "Ter_Mammals"
top_landuse3 <- landuse_change
top_landuse3$taxa <- "Ter_Birds"
top_landuse <- rbind(top_landuse1, top_landuse2, top_landuse3); rm(landuse_change)

# Combine top climate and top landuse
top_landuse$year <- as.numeric(top_landuse$year)
Comb <- dplyr::full_join(climate,top_landuse)

## Get species data
all_taxa_dist_smallrange <- rbind(read.csv("data/amphibians_dist_smallrange.csv"), 
                                  rbind(read.csv("data/birds_dist_smallrange.csv"), 
                                        read.csv("data/mammals_dist_smallrange.csv")))
colnames(all_taxa_dist_smallrange)[5] <- "taxa"
levels(all_taxa_dist_smallrange$taxa) <- c("Amphibians", "Ter_Birds", "Ter_Mammals")
str(all_taxa_dist_smallrange)

# Re-run species data for all models
all_taxa_dist_smallrange <- lapply(c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5"), function(model){
  data <- all_taxa_dist_smallrange
  data$model <- model
  return(data)
})
library(magrittr)
all_taxa_dist_smallrange %<>% bind_rows

# Join taxa and climate/lu but, only keep matching values
smallrange_clim_lu_2050 <- inner_join(all_taxa_dist_smallrange, 
                                            Comb[Comb$year == 2050,])
smallrange_clim_lu_2080 <- inner_join(all_taxa_dist_smallrange, Comb[Comb$year == 2080,])

summary(smallrange_clim_lu_2080)
smallrange_clim_lu <- rbind(smallrange_clim_lu_2050, smallrange_clim_lu_2080)

# Count number of species for each category
smallrange_clim_lu <- smallrange_clim_lu %>% group_by(taxa, species, scenario, year, model) %>% 
  summarise_at(c("dist", "biofuel_cropland", "cropland", "pastures"), sum, na.rm=T) %>%
  mutate_at(c("dist", "biofuel_cropland", "cropland", "pastures"), funs(as.numeric(. > 0)))

# Define categories
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1] <- 1
smallrange_clim_lu$threat[smallrange_clim_lu$biofuel_cropland == 1] <- 2
smallrange_clim_lu$threat[smallrange_clim_lu$cropland == 1] <- 3
smallrange_clim_lu$threat[smallrange_clim_lu$pastures == 1] <- 4
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$biofuel_cropland == 1] <- 5
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$cropland == 1] <- 6
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$pastures == 1] <- 7
smallrange_clim_lu$threat[smallrange_clim_lu$biofuel_cropland == 1 & smallrange_clim_lu$cropland == 1] <- 8
smallrange_clim_lu$threat[smallrange_clim_lu$biofuel_cropland == 1 & smallrange_clim_lu$pastures == 1] <- 9
smallrange_clim_lu$threat[smallrange_clim_lu$cropland == 1 & smallrange_clim_lu$pastures == 1] <- 10
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$biofuel_cropland == 1 & 
                          smallrange_clim_lu$cropland == 1] <- 11
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$biofuel_cropland == 1 & 
                            smallrange_clim_lu$pastures == 1] <- 12
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$cropland == 1 & 
                            smallrange_clim_lu$pastures == 1] <- 13
smallrange_clim_lu$threat[smallrange_clim_lu$biofuel_cropland == 1 & smallrange_clim_lu$cropland == 1 & 
                            smallrange_clim_lu$pastures == 1] <- 14
smallrange_clim_lu$threat[smallrange_clim_lu$dist == 1 & smallrange_clim_lu$biofuel_cropland == 1 & 
                            smallrange_clim_lu$cropland == 1 & smallrange_clim_lu$pastures == 1] <- 15
summary(smallrange_clim_lu$threat)

clim_lu_sum <- smallrange_clim_lu %>% 
  dplyr::select(-dist, -biofuel_cropland, -cropland, -pastures) %>% 
  ungroup() %>% group_by(taxa, scenario, year, threat, model) %>% dplyr::summarise(n = n())

# Change factor labels
clim_lu_sum$taxa <- factor(clim_lu_sum$taxa, labels=c("Amphibians", "Birds", "Mammals"))
clim_lu_sum$scenario <- factor(clim_lu_sum$scenario, labels=c("RCP2.6", "RCP6.0"))

# Re-order var factor
data.frame(clim_lu_sum)

#-#-# Barplot small range species

# Add labels to dataframe
total_clim_lu <- clim_lu_sum %>% group_by(scenario, year, threat, model) %>%
  dplyr::summarise(label = sum(n, na.rm=T))

# Identify number of species under BC and combination with BC
total_clim_lu$threat <- factor(total_clim_lu$threat,
                             labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
                                      "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
                                      "CC+BC+PA", "CC+CR+PA", "BC+CR+PA", "CC+BC+CR+PA"))

total_clim_lu %>% ungroup() %>% filter(threat %in% c("CC")) %>%
  group_by(scenario, year, model) %>% dplyr::summarise(n=sum(label)) %>% 
  group_by(scenario, year) %>%
  dplyr::summarise(mean=mean(n), sd=sd(n))

total_clim_lu %>%   group_by(scenario, year, model) %>% 
  dplyr::summarise(n=sum(label)) %>% 
  group_by(scenario, year) %>%
  dplyr::summarise(mean=mean(n), sd=sd(n))

total_clim_lu %>% ungroup() %>% filter(threat %in% c("BC", "CC+BC", "BC+CR", "BC+PA", "CC+BC+CR", "CC+BC+PA",
                                     "BC+CR+PA", "CC+BC+CR+PA")) %>%
  group_by(scenario, year, model) %>% dplyr::summarise(n=sum(label)) %>% group_by(scenario, year) %>%
  dplyr::summarise(mean=mean(n), sd=sd(n))

# Adapt data for plotting
clim_lu_sum <- left_join(clim_lu_sum, total_clim_lu)
df <- expand.grid(taxa=c("Amphibians", "Birds", "Mammals"),
            scenario = c("RCP2.6", "RCP6.0"), year=c(2050,2080),
            threat=c(1:15))
clim_lu_sum <- left_join(df, clim_lu_sum)
clim_lu_sum$threat <- factor(clim_lu_sum$threat,
                         labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
                                  "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
                                  "CC+BC+PA", "CC+CR+PA", "BC+CR+PA", "CC+BC+CR+PA"))
clim_lu_sum$threat <- factor(clim_lu_sum$threat, levels=levels(clim_lu_sum$threat)[rev(c(1:15))])

# Filter and calculate mean and sd
clim_lu_sum %<>% filter(year == 2080) %>% 
  group_by(taxa, scenario, threat) %>% 
  summarise(mean=mean(n), sd=sd(n))

clim_lu_sum %<>% arrange(desc(taxa))

# Calc error bar position
library(plyr)
plydat <- ddply(clim_lu_sum,.(scenario,threat),transform,
                ystart = cumsum(mean),yend = cumsum(mean) + sd)

# Plot for 2080
ggplot(data=plydat, aes(as.factor(threat), mean)) +   
  geom_bar(aes(fill = taxa), stat="identity") +
  geom_errorbar(aes(ymin=ystart-sd, ymax=yend), width=0.4) + 
  geom_text(aes(y=label, label=label, hjust=-0.5)) +
  facet_wrap(~scenario, nrow=2, scales="free") +
  scale_fill_grey(name="Taxa") + 
  scale_y_continuous(name="Number of species", limits=c(0,2100), expand=c(0,0)) + 
  scale_x_discrete(name="Threat") + coord_flip() + 
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=8), 
        strip.text.x = element_text(size=10, face="bold"),
        strip.background= element_blank(),
        legend.position = c(0.9, 0.13))
ggsave("figures/Figure_4.png", width=6, height=6, dpi=1200, bg="transparent")
