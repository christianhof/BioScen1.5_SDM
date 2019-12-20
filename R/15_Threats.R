# Spatial distribution of climate and land-use threats

# Land-use threat (25% Quantile) was calculated separately

# This code requires the summedProbability, as well as the Landuse Change, files

Colors <- c('red', 'yellow', 'darkgreen', 'darkgoldenrod', 'khaki3', 'pink',
            'limegreen', 'blue','orange','cornflowerblue',
            'black', 'purple','gray60','aquamarine2','brown')

# Extract cells with 25% most change in climatic suitability
# Plot threat maps and calculate total area of each individual threat and threat overlap
library(plyr); library(dplyr); library(ggplot2); library(patchwork)
library(magrittr); library(colorfulVennPlot); library(grid)
threat_areas <- lapply(c("nodisp", "disp"),function(disp){
  year <- 2080
  
  # Load summed Probability 
  sumProb <- lapply(1:3, function(i){
    taxa <- c("Amphibian", "Ter_Bird", "Ter_Mammal")
    sumProb1 <- read.csv(paste0("data/", taxa[i], "_prob_GAM_", disp, ".csv.xz"))
    sumProb2 <- read.csv(paste0("data/", taxa[i], "_prob_GBM_", disp ,".csv.xz"))
    sumProb <- left_join(sumProb1, sumProb2, by=c("x","y"))
    sumProb$taxa <-  c("Amphibians", "Birds", "Mammals")[i]
    return(sumProb)
  })
  sumProb %<>% bind_rows
  
  delta_mean <- sumProb %>% mutate_at(vars("GFDL.ESM2M_rcp26_2050.x":"MIROC5_rcp60_2080.y"), funs(. - EWEMBI_1995.x)) %>% 
    dplyr::select(-c(EWEMBI_1995.x, EWEMBI_1995.y))
  
  delta_rcp26 <- delta_mean %>% select(x,y,taxa,matches(paste0("rcp26_", year))) %>% tidyr::gather(model, value, -c(x,y,taxa))
  delta_rcp26$scenario <- "rcp26"
  delta_rcp60 <- delta_mean %>% select(x,y,taxa,matches(paste0("rcp60_", year))) %>% tidyr::gather(model, value, -c(x,y,taxa))
  delta_rcp60$scenario <- "rcp60"
  delta_all <- bind_rows(delta_rcp26, delta_rcp60); rm(delta_rcp26, delta_rcp60)
  delta_all$year <- year
  delta_all$algorithm <- sapply(delta_all$model, function(x){
    sub("y", "GBM", sub("x", "GAM", strsplit(strsplit(x, split="_")[[1]][3], split="[.]")[[1]][2]))})
  delta_all$model <- sapply(delta_all$model, function(x) strsplit(x, split="_")[[1]][1])
  
  # Calculate ensemble mean
  mean_delta_top25 <- lapply(c("rcp26", "rcp60"), function(x){
    sumProb_sub <- delta_mean %>% select(c(x,y,taxa), matches(paste0(x, "_", year)))
    sumProb_sub$mean <- sumProb_sub %>% select(-c(x,y,taxa)) %>% apply(1, mean, na.rm=TRUE)
    sumProb_sub %<>% select(c(x,y,taxa,mean))
    sumProb_sub$scenario <- x
    return(sumProb_sub)
  })
  mean_delta_top25 %<>% bind_rows
  head(mean_delta_top25)
  
  # Define overall climate change threshhold
  (thr_cc <- delta_all %>% group_by(taxa) %>% 
      filter(value < 0) %>% summarise(thr = quantile(value, probs=0.25)))
  
  # Define climate change threshold for each GCM and algorithm
  #(thr_cc <- delta_all %>% group_by(taxa, model, algorithm) %>% 
  #    filter(value < 0) %>% summarise(thr = quantile(value, probs=0.25)))
  
  # Define mean climate change threshold
  #(thr_cc <- mean_delta_top25 %>% group_by(taxa) %>% 
  #    filter(mean < 0) %>% summarise(thr = quantile(mean, probs=0.25)))
  #ACHTUNG: FUER NEGATIVE WERTE 25%-QUANTILE NEHMEN!
  
  # Add threshhold to data
  delta_all %<>% left_join(thr_cc) %>% filter(value < thr) %>% 
    transmute(x,y,taxa,model,algorithm,scenario,climate = 1)
  mean_delta_top25 %<>% left_join(thr_cc) %>% filter(mean < thr) %>% 
    transmute(x,y,taxa,scenario,climate = 1)
  
  # Read land_use data
  landuse_all <- read.csv("data/landuse_all.csv.xz")
  
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
  top25_lu <- left_join(delta_landuse, threshold) %>% group_by(x,y,scenario, model, var, year) %>% 
    filter(value > threshold) %>% select(-threshold) %>% mutate(value = 1)
  
  # Subset mean data by threshhold
  mean_top25_lu <- left_join(delta_mean, threshold) %>% group_by(x,y,scenario, var, year) %>% 
    filter(mean > threshold) %>% select(-threshold) %>% mutate(mean = 1)
  
  # Delta top25 by GCM and Algorithm
  delta_lu_top25 <- bind_rows(lapply(c("Amphibians", "Birds", "Mammals"), function(taxa){
    data <- top25_lu
    data$taxa <- taxa
    return(data)
  }))
  delta_lu_top25$algorithm <- "GAM"
  delta_lu_top252 <- delta_lu_top25
  delta_lu_top252$algorithm <- "GBM"
  delta_lu_top25 <- bind_rows(delta_lu_top25, delta_lu_top252)
  delta_lu_top25 <- delta_lu_top25[delta_lu_top25$year == year,]
  delta_lu_top25 %<>% ungroup %>% select(-year) %>% tidyr::spread(var, value) 
  delta_all$model <- gsub("[.]", "-", delta_all$model)
  delta_top25 <- full_join(delta_all, delta_lu_top25, by=c("x","y","taxa","model","scenario","algorithm"))
  delta_top25$biofuel_cropland[is.na(delta_top25$biofuel_cropland)] <- 0
  delta_top25$cropland[is.na(delta_top25$cropland)] <- 0
  delta_top25$pastures[is.na(delta_top25$pastures)] <- 0
  delta_top25$climate[is.na(delta_top25$climate)] <- 0
  delta_top25 <- tidyr::gather(delta_top25, var, value, -c(x,y,taxa,scenario,model,algorithm))
  delta_top25$var <- factor(delta_top25$var, labels=c("Biofuel cropland", "Climate", 
                                                      "Non-biofuel cropland", "Pastures"))
  delta_top25$taxa <- factor(delta_top25$taxa, labels=c("Amphibians", "Birds", "Mammals"))
  delta_top25$scenario <- factor(delta_top25$scenario, labels=c("RCP2.6", "RCP6.0"))
  
  # Mean delta top25
  delta_lu_top25 <- bind_rows(lapply(c("Amphibians", "Birds", "Mammals"), function(taxa){
    data <- mean_top25_lu
    data$taxa <- taxa
    return(data)
  }))
  delta_lu_top25 <- delta_lu_top25[delta_lu_top25$year == year,]
  delta_lu_top25 %<>% ungroup() %>% select(-year) %>% tidyr::spread(var, mean)
  mean_top25 <- full_join(mean_delta_top25, delta_lu_top25)
  mean_top25$biofuel_cropland[is.na(mean_top25$biofuel_cropland)] <- 0
  mean_top25$cropland[is.na(mean_top25$cropland)] <- 0
  mean_top25$pastures[is.na(mean_top25$pastures)] <- 0
  mean_top25$climate[is.na(mean_top25$climate)] <- 0
  mean_top25 <- tidyr::gather(mean_top25, var, value, -c(x,y,taxa,scenario))
  mean_top25$var <- factor(mean_top25$var, labels=c("Biofuel cropland", "Climate", 
                                                    "Non-biofuel cropland", "Pastures"))
  mean_top25$taxa <- factor(mean_top25$taxa, labels=c("Amphibians", "Birds", "Mammals"))
  mean_top25$scenario <- factor(mean_top25$scenario, labels=c("RCP2.6", "RCP6.0"))
  
  # Remove 0 values from data frame
  delta_top25$value[delta_top25$value == 0] <- NA
  mean_top25$value[mean_top25$value == 0] <- NA
  
  # Turn data into wide format
  delta_top25  %<>% tidyr::spread(var, value)
  mean_top25 %<>% tidyr::spread(var, value)
  
  #Remove rows where all variables are NA
  delta_top25 <- delta_top25[rowSums(is.na(delta_top25)) < 4,]
  mean_top25 <- mean_top25[rowSums(is.na(mean_top25)) < 4,]
  
  # Define categories
  delta_top25$one[delta_top25$`Climate` == 1] <- 1
  delta_top25$two[delta_top25$`Biofuel cropland` == 1] <- 2
  delta_top25$three[delta_top25$`Non-biofuel cropland` == 1] <- 3
  delta_top25$four[delta_top25$`Pastures` == 1] <- 4
  
  delta_top25$twelve[delta_top25$`Climate` == 1 & delta_top25$`Biofuel cropland` == 1] <- 12
  delta_top25$thirteen[delta_top25$`Climate` == 1 & delta_top25$`Non-biofuel cropland` == 1] <- 13 
  delta_top25$fourteen[delta_top25$`Climate` == 1 & delta_top25$`Pastures` == 1] <- 14
  delta_top25$twentythree[delta_top25$`Biofuel cropland` == 1 & delta_top25$`Non-biofuel cropland` == 1] <- 23
  delta_top25$twentyfour[delta_top25$`Biofuel cropland` == 1 & delta_top25$`Pastures` == 1] <- 24
  delta_top25$thirtyfour[delta_top25$`Non-biofuel cropland` == 1 & delta_top25$`Pastures` == 1] <- 34
  
  delta_top25$onetwothree[delta_top25$`Climate` == 1 & delta_top25$`Biofuel cropland` == 1 & 
                            delta_top25$`Non-biofuel cropland` == 1] <- 123
  delta_top25$onetwofour[delta_top25$`Climate` == 1 & delta_top25$`Biofuel cropland` == 1 & 
                           delta_top25$`Pastures` == 1] <- 124
  delta_top25$onethreefour[delta_top25$`Climate` == 1 & delta_top25$`Non-biofuel cropland` == 1 & 
                             delta_top25$`Pastures` == 1] <- 134
  delta_top25$twothreefour[delta_top25$`Biofuel cropland` == 1 & delta_top25$`Non-biofuel cropland` == 1 & 
                             delta_top25$`Pastures` == 1] <- 234
  delta_top25$onetwothreefour[rowSums(is.na(delta_top25)) == 0] <- 1234
  delta_top25[, "max"] <- delta_top25 %>% select(one:onetwothreefour) %>% apply(1, max, na.rm=T)
  
  # Calculate dominant threat
  #mean_top25 <- delta_top25 %>% select(x,y,taxa,model,algorithm,scenario,max) %>% 
  #  group_by(x,y,taxa,scenario,max) %>%
  #  summarise (n = n()) %>%
  #  mutate(freq = (n/8)*100) %>% filter(freq > 50)
  
  # Define categories of ensemble mean threat
  mean_top25$one[mean_top25$`Climate` == 1] <- 1
  mean_top25$two[mean_top25$`Biofuel cropland` == 1] <- 2
  mean_top25$three[mean_top25$`Non-biofuel cropland` == 1] <- 3
  mean_top25$four[mean_top25$`Pastures` == 1] <- 4
  
  mean_top25$twelve[mean_top25$`Climate` == 1 & mean_top25$`Biofuel cropland` == 1] <- 12
  mean_top25$thirteen[mean_top25$`Climate` == 1 & mean_top25$`Non-biofuel cropland` == 1] <- 13 
  mean_top25$fourteen[mean_top25$`Climate` == 1 & mean_top25$`Pastures` == 1] <- 14
  mean_top25$twentythree[mean_top25$`Biofuel cropland` == 1 & mean_top25$`Non-biofuel cropland` == 1] <- 23
  mean_top25$twentyfour[mean_top25$`Biofuel cropland` == 1 & mean_top25$`Pastures` == 1] <- 24
  mean_top25$thirtyfour[mean_top25$`Non-biofuel cropland` == 1 & mean_top25$`Pastures` == 1] <- 34
  
  mean_top25$onetwothree[mean_top25$`Climate` == 1 & mean_top25$`Biofuel cropland` == 1 & 
                           mean_top25$`Non-biofuel cropland` == 1] <- 123
  mean_top25$onetwofour[mean_top25$`Climate` == 1 & mean_top25$`Biofuel cropland` == 1 & 
                          mean_top25$`Pastures` == 1] <- 124
  mean_top25$onethreefour[mean_top25$`Climate` == 1 & mean_top25$`Non-biofuel cropland` == 1 & 
                            mean_top25$`Pastures` == 1] <- 134
  mean_top25$twothreefour[mean_top25$`Biofuel cropland` == 1 & mean_top25$`Non-biofuel cropland` == 1 & 
                            mean_top25$`Pastures` == 1] <- 234
  mean_top25$onetwothreefour[rowSums(is.na(mean_top25)) == 0] <- 1234
  mean_top25[, "max"] <- mean_top25 %>% select(one:onetwothreefour) %>% apply(1, max, na.rm=T)
  
  rm(delta_all, delta_landuse, delta_lu_top25, delta_lu_top252, delta_mean, landuse_all,
     mean_delta_top25, mean_top25_lu, sumProb, thr_cc, threshold, top25_lu)
  
  ###########################
  
  # Create plot
  data(outline, package="ggmap2")
  outline <- sf::st_as_sf(outline)
  
  # Load animals and plot them!
  bird <- raster::stack("data/bird.png")
  mammal <- raster::stack("data/mammal.png")
  toad <- raster::stack("data/toad.png")
  
  library(RStoolbox)
  p3 <- ggRGB(toad) + theme_classic() + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),
          axis.text = element_blank(), axis.line = element_blank())
  p4 <- ggRGB(bird) + theme_classic() + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),
          axis.text = element_blank(), axis.line = element_blank())
  p5 <- ggRGB(mammal) + theme_classic() + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title = element_blank(), axis.ticks=element_blank(),
          axis.text = element_blank(), axis.line = element_blank())
  
  xlim <- c(floor(min(mean_top25$x)), ceiling(max(mean_top25$x)))
  ylim <- c(floor(min(mean_top25$y)), ceiling(max(mean_top25$y)))
  mean_top25$max <- factor(mean_top25$max, 
                           labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
                                    "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
                                    "CC+BC+PA", "CC+CR+PA", "BC+CR+PA", "CC+BC+CR+PA"))
  #Dominant threat has less categories
  #mean_top25$max <- factor(mean_top25$max, 
  #                         labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
  #                                  "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
  #                                  "CC+BC+PA","BC+CR+PA"))
  
  p1 <- ggplot() + geom_tile(data=mean_top25, aes(x=x, y=y, fill=max)) + 
    facet_grid(taxa ~ scenario) + 
    geom_sf(data=outline, fill="transparent", colour="black", size=0.25) + 
    scale_fill_manual(name="", na.value="transparent", values=Colors) + 
    theme_classic() + theme(axis.line = element_blank(), axis.ticks = element_blank(), 
                            axis.text = element_blank(), axis.title = element_blank(),
                            strip.text.x = element_text(size = 16, face= "bold"),
                            strip.text.y = element_blank(), strip.placement = "outside", 
                            strip.background= element_blank(), legend.position="none") + 
    coord_sf(expand=F, xlim=c(-180,180), ylim=c(-65,85), ndiscr=0)
  
  # Calculate area of each category
  data(landseamask_generic, package="rISIMIP")
  area <- raster::area(landseamask_generic, na.rm=TRUE)
  area <- data.frame(raster::rasterToPoints(area))
  colnames(area) <- c("x", "y", "area")
  delta_top25 %<>% left_join(area)
  delta_top25 %<>% dplyr::select(scenario, taxa, area, model, algorithm, max) %>% 
    group_by(scenario, taxa, model, algorithm, max) %>%
    summarise(sum=sum(area, na.rm=T)) %>% ungroup()
  delta_top25$dispersal <- disp
  delta_top25$scenario <- factor(delta_top25$scenario, labels=c("RCP2.6", "RCP6.0"))
  delta_top25 <- tidyr::drop_na(delta_top25)
  colnames(delta_top25) <- c("Scenario", "Taxa", "Model", "Algorithm", "Threat", "Area_km2", "Dispersal")
  delta_top25$area_perc <- round(delta_top25$`Area_km2`/sum(area$area)*100,digits=2)
  delta_top25 <- delta_top25 %>% mutate(Area_Mha = round(Area_km2/10000,digits=2))
  
  # Calculate mean and error and create plot
  delta_sum <- delta_top25 %>% group_by(Scenario, Taxa, Threat) %>% 
    dplyr::summarise(mean=mean(area_perc), sd=sd(area_perc))
  delta_sum$Threat <- factor(delta_sum$Threat, 
                             labels=c("CC", "BC", "CR", "PA", "CC+BC", "CC+CR", 
                                      "CC+PA", "BC+CR", "BC+PA", "CR+PA", "CC+BC+CR", 
                                      "CC+BC+PA", "CC+CR+PA", "BC+CR+PA", "CC+BC+CR+PA"))
  delta_sum %<>% arrange(desc(Threat))
  
  # Calc Error bar position
  plydat <- ddply(delta_sum,.(Taxa,Scenario),transform,
                  ystart = cumsum(mean),yend = cumsum(mean) + sd)
  
  p2 <- ggplot(data=plydat, aes(x=Taxa, y=mean, fill=Threat)) + 
    geom_bar(stat="identity", width=0.55, color=NA) + 
    facet_wrap(~ Scenario, scales="free", ncol=1) + 
    #geom_segment(aes(x = Taxa, xend=Taxa, y = ystart,yend = yend)) + 
    #geom_point(aes(x = Taxa,y = yend), shape = "|",show.legend = FALSE) + coord_flip() + 
    geom_errorbar(aes(ymin=ystart-sd, ymax=yend), width=0.1) + 
    scale_fill_manual(values=c('red', 'yellow', 'darkgreen', 'darkgoldenrod', 'khaki3', 'pink',
                               'limegreen', 'blue','orange','cornflowerblue',
                               'black', 'purple','gray60','aquamarine2','brown')) + 
    scale_y_continuous(limits=c(0,45), expand=c(0,0)) + 
    labs(x="", y="% of land surface") + theme_classic() + 
    theme(axis.title = element_text(size = 14, face= "bold"), legend.position="none",
          axis.text.x = element_blank(), axis.text.y = element_text(size=12),
          strip.text = element_text(size = 14, face= "bold"),
          strip.background= element_blank())
  p2 <- {p2 + plot_spacer() + plot_layout(ncol=1, heights=c(2,1)) + theme_classic()}
  
  # Infos for Venn diagram
  regions <- rep("", 15) #regions <- seq(15) #If you want to and numbers to each colour
  names(regions) <- c('1000', '0100', '0010', '0001', '1100', '1010', 
                      '1001', '0110','0101', '0011', 
                      '1110', '1101', '1011', '0111', '1111')
  vp1 <- viewport(width=0.7, height=1, x = 0, y = 0, just=c("left","bottom"))
  vp2 <- viewport(width=0.13, height=0.13, x = -0.03, y = 0.65, just=c("left", "bottom"))
  vp3 <- viewport(width=0.13, height=0.13, x = -0.03, y = 0.35, just=c("left", "bottom"))
  vp4 <- viewport(width=0.13, height=0.13, x = -0.03, y = 0.07, just=c("left", "bottom"))
  vp5 <- viewport(width=0.15, height=0.28, x = 0.8, y = 0, just=c("left", "bottom"))
  vp6 <- viewport(width=0.1, height=0.1, x = 0.74, y = 0.25, just=c("left", "bottom"))
  vp7 <- viewport(width=0.1, height=0.1, x = 0.815, y = 0.25, just=c("left", "bottom"))
  vp8 <- viewport(width=0.1, height=0.1, x = 0.895, y = 0.25, just=c("left", "bottom"))
  vp9 <- viewport(width=0.3, height=0.99, x = 0.7, y = 0, just=c("left", "bottom"))
  
  if(disp == "nodisp" & year == 2080){
    # Create map with Venn diagram
    cairo_pdf(file=paste0("figures/Figure_2.pdf"), 
              width=12, height=6, bg="transparent", fallback_resolution = 900)
    plot.new()
    print(p1, vp=vp1)
    print(p2, vp=vp9)
    print(p3, vp=vp2)
    print(p4, vp=vp3)
    print(p5, vp=vp4)
    print(p3, vp=vp6)
    print(p4, vp=vp7)
    print(p5, vp=vp8)
    pushViewport(vp5)
    plotVenn4d(regions, Colors=Colors, labels=c("CC", "BC", "CR", "PA"), Title = '')
    dev.off()
  } else{
    png(file=paste0("figures/threat_map_disp.png"), 
              width=12, height=6, bg="transparent", units="in", res = 900)
    plot.new()
    print(p, vp=vp1)
    print(p3, vp=vp2)
    print(p4, vp=vp3)
    print(p5, vp=vp4)
    print(p3, vp=vp6)
    print(p4, vp=vp7)
    print(p5, vp=vp8)
    pushViewport(vp5)
    plotVenn4d(regions, Colors=Colors, labels=c("CC", "BC", "CR", "PA"), Title = '')
    dev.off()
  }
  delta_top25$year <- year
  delta_top25$disp <- disp
  return(delta_top25)
})
threat_areas %<>% bind_rows
colnames(threat_areas)
library(tidyr)

(threat_areas_Mha <- threat_areas %>% select(-area_perc, -Area_km2, -disp) %>% 
  unite("disp_scenario", Dispersal, Scenario) %>% 
  group_by(disp_scenario, Taxa, Threat, year) %>% 
  dplyr::summarise(mean=round(mean(Area_Mha),2), sd=round(sd(Area_Mha),2)) %>% 
  unite(mean, sd, col="Area_Mha", sep=" ± ") %>%
  spread(disp_scenario, Area_Mha))
colnames(threat_areas_Mha)[4:7] <- c("disp_RCP2.6_area_Mha", "disp_RCP6.0_area_Mha", 
                                     "nodisp_RCP2.6_area_Mha", "nodisp_RCP6.0_area_Mha")
(threat_areas_perc <- threat_areas %>% select(-Area_Mha, -Area_km2, -disp) %>% 
  unite("disp_scenario", Dispersal, Scenario) %>% 
  group_by(disp_scenario, Taxa, Threat, year) %>% 
  dplyr::summarise(mean=round(mean(area_perc),2), sd=round(sd(area_perc),2)) %>% 
  unite(mean, sd, col="area_perc", sep=" ± ") %>%
  spread(disp_scenario, area_perc))
colnames(threat_areas_perc)[4:7] <- c("disp_RCP2.6_area_perc", "disp_RCP6.0_area_perc", 
                                      "nodisp_RCP2.6_area_perc", "nodisp_RCP6.0_area_perc")
threat_areas_all <- full_join(threat_areas_Mha, threat_areas_perc)
colnames(threat_areas_all)
threat_areas_all <- threat_areas_all[,c(1,2,3,6,10,7,11,4,8,5,9)] %>% 
  arrange(year, Taxa, `Threat`)
threat_areas_2080 <- threat_areas_all %>% filter(year == 2080) %>% select(-year)
write.csv2(threat_areas_2080, paste0("data/threat_areas_2080.csv"), row.names=F)

threat_areas_2080 <- read_delim("data/threat_areas_2080.csv", 
                                    +     ";", escape_double = FALSE, trim_ws = TRUE)

threat_areas_2080 %>% separate(nodisp_RCP2.6_area_perc, c("perc", "sd"), sep=" ± ") %>%  
  mutate(perc=as.numeric(perc), sd=as.numeric(sd)) %>% 
  group_by(Taxa) %>% summarise(perc=sum(perc), sd=sum(sd))
threat_areas_2080 %>% separate(nodisp_RCP6.0_area_perc, c("perc", "sd"), sep=" ± ") %>%  
  mutate(perc=as.numeric(perc), sd=as.numeric(sd)) %>% 
  group_by(Taxa) %>% summarise(perc=sum(perc), sd=sum(sd))

threat_areas_2080 %>% 
  separate(nodisp_RCP2.6_area_perc, c("perc", "sd"), sep=" ± ") %>%  
  mutate(perc=as.numeric(perc), sd=as.numeric(sd)) %>% 
  filter(Threat %in% c(12,13,123,124,134,1234)) %>% 
  group_by(Taxa) %>% summarise(perc=sum(perc), sd=sum(sd))
threat_areas_2080 %>% 
  separate(nodisp_RCP6.0_area_perc, c("perc", "sd"), sep=" ± ") %>%  
  mutate(perc=as.numeric(perc), sd=as.numeric(sd)) %>% 
  filter(Threat %in% c(12,13,123,124,134,1234)) %>% 
  group_by(Taxa) %>% summarise(perc=sum(perc), sd=sum(sd))

