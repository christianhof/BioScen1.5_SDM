# Create prediction maps for the different taxa, GCMs, rcps and timesteps

# Empty R environment
rm(list=ls())

# Load libraries
packages <- c("dplyr", "snowfall", "data.table")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

# Set file directory
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/bigdata/mbiber/data" # darkstar
filedir <- "/bigdata_local/mbiber" # ceremony - mammals
#filedir <- "/home/mbiber/data" # ceremony - birds

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 2

# Model type
k <- 4; model_type <- c("GAM", "GBM", "MaxEnt", "RF")[k]

# Final species names
spNames <- sapply(list.files(paste0(filedir, "/", taxa[i], "_", 
                                    model_type, "_results_climate")), function(x){
                                      paste0(strsplit(x, split="_")[[1]][1:2], collapse="_")
                                    })

# Final prediction files
predFiles <- list.files(paste0(filedir, "/", taxa[i], "_", model_type, "_results_climate"), 
                        full.names=T)
predFiles <- sapply(spNames, function(x) predFiles[grepl(predFiles, pattern=x)][[1]])
length(predFiles)

# Read csv files
sfInit(parallel=TRUE, cpus=ceiling(0.75*parallel::detectCores()))
predData <- sfLapply(predFiles, function(x){readr::read_csv(x)})
sfStop()
#predData <- dplyr::bind_rows(predData)
predData <- data.table::rbindlist(predData)

# Calculate summed probability and save to .csv.xz file
lapply(c("presence", "dispersal1", "dispersal2", "dispersal3", "dispersal4", "fulldisp"), function(disp){
  if(!file.exists(paste0(filedir, "/", taxa[i], "_prob_", 
                         model_type, "_", disp, ".csv.xz"))){
    if(disp != "fulldisp"){
      sumProb_nodisp <- predData %>% group_by(x, y) %>% filter(1 == !!as.name(disp)) %>% 
        select(x,y,GFDL.ESM2M_piControl_1845:MIROC5_rcp26_2250) %>% 
        mutate_if(is.character, as.numeric) %>% 
        summarise_all(sum, na.rm=T) %>% as.data.frame()
      readr::write_csv(sumProb_nodisp, paste0(filedir, "/", taxa[i], "_prob_", 
                                              model_type, "_", disp, ".csv.xz"))
    } else{
      sumProb_disp <- predData %>% group_by(x, y) %>% 
        select(x,y,GFDL.ESM2M_piControl_1845:MIROC5_rcp26_2250) %>% 
        mutate_if(is.character, as.numeric) %>% 
        summarise_all(sum, na.rm=T) %>% as.data.frame()
      readr::write_csv(sumProb_disp, paste0(filedir, "/", taxa[i], "_prob_", 
                                            model_type, "_fulldisp.csv.xz"))
    }
  }
  return(NULL)
}); q(save="no")

##############################

# Save summed probability to .nc file

# Empty R environment
rm(list=ls())

# Load libraries
packages <- c("dplyr", "snowfall", "tidyr", "ncdf4", "readr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

# Set file directory
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/scratch/mbiber/data" # darkstar
#filedir <- "/bigdata_local/mbiber" # ceremony
filedir <- "/home/matt/Documents"

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 1

# Model type
k <- 1; model_type <- c("GAM", "GBM", "MaxEnt", "RF")[k]

# Set up NetCDF file combinations
df <- expand.grid(rcp=c("rcp26", "rcp60", "historical", "piControl"), 
                  model=c("MIROC5", "HadGEM2.ES", "IPSL.CM5A.LR", "GFDL.ESM2M"))
df <- rbind(expand.grid(rcp="1995", model="EWEMBI"), df)
df <- rbind(df, expand.grid(rcp="2100rcp26", 
                            model=c("MIROC5", "HadGEM2.ES", "IPSL.CM5A.LR")))

library(snowfall)
sfInit(parallel=TRUE, cpus=10)
sfLibrary(dplyr); sfLibrary(tidyr); sfLibrary(ncdf4)
sfExport(list=c("filedir", "taxa", "model_type", "i", "df"))
sfLapply(1:nrow(df), function(rcp_mod){
  #Define years
  if(df$rcp[rcp_mod] == "1995"){
    year <- 1995
  } else if(df$rcp[rcp_mod] == "historical"){
    year <- 1990
  } else if(df$rcp[rcp_mod] == "piControl"){
    year <- 1845
  } else if(df$rcp[rcp_mod] == "rcp26"){
    year <- c(2009, 2010, 2020, 2026, 2032, 2048, 2050, 2052, 2056, 2080)
  } else if(df$rcp[rcp_mod] == "rcp60"){
    year <- c(2009, 2010, 2020, 2026, 2032, 2048, 2050, 2052, 2056, 2080)
  } else if(df$rcp[rcp_mod] == "2100rcp26"){
    year <- c(2100, 2150, 2200, 2250)
  }
  if(any(year == 1995)){
    filename <- paste0(filedir, "/OutputData/biodiversity/bioscen1.5-sdm-", 
                       tolower(model_type),
                       "_ewembi_nobc_hist_nosoc_co2_", tolower(taxa[i]), 
                       "sumprob_global_30year-mean_", 
                       min(year), "_", max(year), ".nc4")
  } else{
    if(df$rcp[rcp_mod] == "2100rcp26"){rcp <- "rcp26"} 
    else{rcp <- df$rcp[rcp_mod]}
    filename <- paste0(filedir, "/OutputData/biodiversity/bioscen1.5-sdm-", 
                       tolower(model_type), "_",
                       tolower(gsub("[.]", "-", df$model[rcp_mod])), 
                       "_ewembi_", rcp, "_nosoc_co2_", 
                       tolower(taxa[i]), "sumprob_global_30year-mean_", 
                       min(year), "_", max(year), ".nc4")
  }
  
  # Define the dimensions
  dimX = ncdim_def(name="lon", units="degrees", vals=seq(-179.75, 179.75, length = 720))
  dimY = ncdim_def(name="lat", units="degrees", vals=seq(89.75, -89.75, length = 360))
  dimT = ncdim_def(name="time", units="years since 1661-1-1 00:00:00", 
                   vals=c(year-1661), calendar="standard")
  
  # Define data for NetCDF file
  vard <- lapply(c("presence", "dispersal1", "dispersal2", 
                   "dispersal3", "dispersal4", "fulldisp"), function(disp){
                     ncvar_def(disp, "Summed probability of occurrence per cell", 
                               list(dimX,dimY,dimT), 1.e+20, prec="double", compression=9)})
  # Create the NetCDF file
  #file.exists(filename)
  nc <- nc_create(filename, vard)
  ncatt_put(nc, varid=0, attname="contact", attval="Matthias Biber <matthias.biber@tum.de>")
  ncatt_put(nc, varid=0, attname="institution", attval="Technical University Munich (Germany)")
  
  lapply(1:6, function(x){
    disp <- c("presence", "dispersal1", "dispersal2", 
              "dispersal3", "dispersal4", "fulldisp")[x]
    
    # Read data
    predData <- readr::read_csv(paste0(filedir, "/", taxa[i], "_prob_", 
                                       model_type, "_", disp, ".csv.xz"))
    
    ## Select all data for one model_rcp combination
    if(df$rcp[rcp_mod] == "2100rcp26"){
      predData <- predData %>% select(x,y, HadGEM2.ES_rcp26_2100:MIROC5_rcp26_2250) %>%
        dplyr::select(x, y, matches(paste(df$model[rcp_mod], "rcp26", sep="_")))
    } else{
      predData <- predData %>% select(x,y, GFDL.ESM2M_piControl_1845:MIROC5_rcp60_2080) %>%
        dplyr::select(x, y, matches(paste(df$model[rcp_mod], 
                                          df$rcp[rcp_mod], sep="_")))
    }
    
    # Expand dataframe with NAs
    df_spat <- expand.grid(x=seq(-179.75, 179.75, length = 720),
                           y=seq(89.75, -89.75, length = 360))
    predData <- left_join(df_spat, predData); rm(df_spat)
    
    # Turn data into array
    predData <- dplyr::select(predData, -x, -y)
    colnames(predData) <- year
    data <- array(unlist(predData), dim=c(720, 360, ncol(predData)), 
                  dimnames=list(NULL, NULL, colnames(predData)))
    
    # Write data to the NetCDF file
    ncvar_put(nc, vard[[x]], data, start=c(1,1,1), count=c(-1,-1,-1))
  })
  # Close your new file to finish writing
  nc_close(nc); rm(nc)
})
sfStop()

# Test NetCDF file
library(raster)
#test <- stack("E://OutputData/biodiversity/
#              bioscen1.5-sdm_EWEMBI_historical_amphibian-sumprob-nodisp.nc4")
#plot(test[[1]])

########################################

# Calculate change with histogram in probability for disp and year
library(dplyr); library(ggplot2); library(patchwork); library(magrittr)
lapply(c("2050", "2080"), function(year){
  time_rcp <- expand.grid(time=year, rcp= c("rcp26", "rcp60")) %>% 
    tidyr::unite("time_rcp", c(rcp, time))
  time_rcp <- as.vector(time_rcp$time_rcp)
  
  delta_mean <- lapply(c("nodisp", "disp"), function(disp){
    
    sumProb <- lapply(1:3, function(i){
      taxa <- c("Amphibian", "Ter_Bird", "Ter_Mammal")
      sumProb1 <- read.csv(paste0("data/", taxa[i], "_prob_GAM_", disp, ".csv.xz"))
      sumProb2 <- read.csv(paste0("data/", taxa[i], "_prob_GBM_", disp ,".csv.xz"))
      sumProb <- left_join(sumProb1, sumProb2, by=c("x","y"))
      sumProb$taxa <-  c("Amphibians", "Birds", "Mammals")[i]
      return(sumProb)
    })
    sumProb %<>% bind_rows 
    
    # Summarise by taxa
    sumtaxa <- sumProb %>% group_by(x,y) %>% 
      dplyr::select(-taxa) %>% summarise_all(sum)
    cv_taxa <- lapply(time_rcp, function(x){
      sumProb_sub <- sumtaxa %>% dplyr::select(c(x,y), matches(x))
      sumProb_sub$CV <- apply(sumProb_sub[,-c(1,2)], 1, raster::cv, na.rm=TRUE)
      sumProb_sub %<>% dplyr::select(c(x,y,CV))
      sumProb_sub$time_rcp <- x
      return(sumProb_sub)
    })
    cv_taxa %<>% bind_rows()
    cv_taxa$time_rcp <- factor(cv_taxa$time_rcp, labels=c(paste0(year, " RCP2.6"), 
                                                          paste0(year, " RCP6.0")))
    cv_taxa$x <- round(cv_taxa$x/2.5)*2.5
    cv_taxa$y <- round(cv_taxa$y/2.5)*2.5
    cv_taxa %<>% group_by(x,y,time_rcp) %>% summarise(CV=max(CV)) %>% filter(CV <= 10)
    #cv_taxa %<>% group_by(x,y,time_rcp) %>% summarise(CV=max(CV)) %>% filter(CV <= 25)
    
    deltaProb <- sumProb %>% 
      mutate_at(vars("GFDL.ESM2M_rcp26_2050.x":"MIROC5_rcp60_2080.y"), 
                funs(. - EWEMBI_1995.x)) %>% 
      dplyr::select(-c(EWEMBI_1995.x, EWEMBI_1995.y))
    deltaChange <- sumtaxa %>% mutate_at(vars("GFDL.ESM2M_rcp26_2050.x":"MIROC5_rcp60_2080.y"), 
                                       funs(. - EWEMBI_1995.x)) %>% 
      dplyr::select(-c(EWEMBI_1995.x, EWEMBI_1995.y))

    # Calculate sign of change
    sign_change <- lapply(time_rcp, function(x){
      deltaProb_sub <- deltaChange %>% dplyr::select(c(x,y), matches(x)) %>% 
        tidyr::gather(var, value, -c(x,y)) %>% filter(value > 0) %>% group_by(x,y) %>% 
        summarise (n = n()) %>% mutate(freq = (n/8)*100)
      deltaProb_sub$time_rcp <- x
     return(deltaProb_sub)
    })
    sign_change %<>% bind_rows
    sign_change$time_rcp <- factor(sign_change$time_rcp, labels=c(paste0(year, " RCP2.6"), 
                                                                  paste0(year, " RCP6.0")))
    
    # Select 5 out of 8 and reduce resolution
    #sign_sub <- sign_change %>% ungroup() %>%
    #  mutate(x=round(x/2.5)*2.5, y=round(y/2.5)*2.5) %>% group_by(x,y,time_rcp) %>% 
    #  summarise(max=max(freq), min=min(freq)) %>% filter(min >= 62.5 | max <= 37.5) 
    
    # Select 6 out of 8 and reduce resolution
    sign_sub <- sign_change %>% ungroup() %>%
      mutate(x=round(x/2.5)*2.5, y=round(y/2.5)*2.5) %>% group_by(x,y,time_rcp) %>% 
      summarise(max=max(freq), min=min(freq)) %>% filter(min >= 75 | max <= 25) 
    
    # Calculate delta mean
    delta_mean <- lapply(time_rcp, function(x){
      deltaProb_sub <- deltaProb %>% dplyr::select(c(x,y,taxa), matches(x))
      deltaProb_sub$mean <- deltaProb_sub %>% dplyr::select(-c("x","y", "taxa")) %>% 
        apply(1, mean, na.rm=TRUE)
      deltaProb_sub %<>% dplyr::select(c(x,y,taxa,mean))
      deltaProb_sub$time_rcp <- x
      return(deltaProb_sub)
    })
    delta_mean %<>% bind_rows
    delta_mean$time_rcp <- factor(delta_mean$time_rcp, labels=c(paste0(year, " RCP2.6"), 
                                                          paste0(year, " RCP6.0")))
    
    delta_all <- delta_mean %>% group_by(x,y,time_rcp) %>% 
      summarise(total=sum(mean))
    delta_x <- delta_mean %>% group_by(x, time_rcp, taxa) %>% 
      summarise(mean=mean(mean, na.rm=T))
    delta_y <- delta_mean %>% group_by(y, time_rcp, taxa) %>% 
      summarise(mean=mean(mean, na.rm=T))
    data(outline, package="ggmap2")
    library(sf)
    outline <- sf::st_as_sf(outline)
    col_val <- scales::rescale(unique(c(seq(min(delta_all$total), 0, length=3), 
                                        seq(0, max(delta_all$total), length=3))))
    lim_map <- c(min(delta_all$total), max(delta_all$total))
    lim_histx <- unlist(c(delta_x %>% group_by(taxa) %>% summarise(min=min(mean)) %>% 
                            ungroup() %>% summarise(sum(min)), 
                          delta_x %>% group_by(taxa) %>% summarise(max=max(mean)) %>% 
                            ungroup() %>% summarise(sum(max))))
    lim_histy <- unlist(c(delta_y %>% group_by(taxa) %>% summarise(min=min(mean)) %>% 
                            ungroup() %>% summarise(sum(min)), 
                          delta_y %>% group_by(taxa) %>% summarise(max=max(mean)) %>% 
                            ungroup() %>% summarise(sum(max))))
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + theme_classic() + 
      theme(legend.position="none", axis.title = element_blank(),
            axis.line = element_blank(), axis.ticks = element_blank(), 
            axis.text = element_blank())
    p2 <- lapply(c("RCP2.6", "RCP6.0"), function(rcp){
      ggplot() + geom_histogram(data=delta_x[delta_x$time_rcp == paste(year, rcp),], 
                                aes(x=x, y=mean, fill=taxa), width=1,
                                stat="identity", position="stack", colour=NA) + 
        scale_y_continuous(position = "right", expand=c(0,0), limits=lim_histx,
                           breaks=c(-60, -30, 0), labels=c(-60, -30, 0)) + 
        scale_x_continuous(expand=c(0,0)) + scale_fill_grey() + theme_classic() + 
        theme(legend.position="none", axis.title = element_blank(),
              axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), axis.text.y = element_text(size=16))
    })
    
    p3 <- lapply(c("RCP2.6", "RCP6.0"), function(rcp){
      ggplot() + geom_histogram(data=delta_y[delta_y$time_rcp == paste(year, rcp),], 
                                aes(x=y, y=mean, fill=taxa), width=1,
                                stat="identity", position="stack", colour=NA) + 
        scale_y_reverse(limits=rev(lim_histy), breaks=c(0,-25, -50), labels=c(0, "", -50)) + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_fill_grey(name="Taxa", labels=c("Amphibians", "Birds", "Mammals")) + 
        coord_flip(expand=FALSE) + theme_classic() + 
        theme(legend.position="bottom", legend.title=element_text(size=16, face="bold"),
              legend.text=element_text(size=16),
              plot.background = element_rect(fill = "transparent"), 
              axis.title = element_blank(), axis.line.y = element_blank(),
              axis.ticks.y = element_blank(), axis.text.y = element_blank(),
              axis.text.x = element_text(size=16)) + 
        guides(fill = guide_legend(direction = "vertical"))
    })

    library(cowplot)
    legend <- ggdraw(get_legend(p3[[1]]))
    
    p4 <- lapply(c("RCP2.6", "RCP6.0"), function(rcp){
      data <- delta_all[delta_all$time_rcp == paste(year, rcp),]
      # Make map
      ggplot() + geom_tile(data=data, aes(x=x, y=y, fill=total)) + 
        #geom_point(data=sign_sub, aes(x=x,y=y), colour="grey50", fill="transparent", size=0.03) + 
        geom_point(data=cv_taxa, aes(x=x,y=y), colour="grey50", fill="transparent", size=0.03) + 
        geom_sf(data=outline, fill="transparent", colour="black") + 
        scale_fill_gradientn(name="", 
                             colours=rev(colorRampPalette(
                               c("#00007F", "blue", "#007FFF", "cyan", 
                                 "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
                             na.value= "grey50", values=col_val, limits=lim_map) + 
        coord_sf(expand=F, xlim=c(min(data$x), max(data$x)), ylim=c(min(data$y),max(data$y)), 
                 ndiscr=0) + theme_classic() + 
        theme(axis.title = element_blank(),axis.line = element_blank(),
              axis.ticks = element_blank(), axis.text = element_blank(),
              plot.background = element_rect(fill = "transparent"), 
              legend.background = element_rect(fill = "transparent"), 
              legend.key.width=unit(2, "cm"), legend.position="bottom",
              legend.text = element_text(size=16), 
              legend.box.background = element_rect(fill = "transparent", colour=NA))
    })
    
    legend2 <- ggdraw(get_legend(p4[[1]]))
    
    pall <- {{{p + ggtitle("a)") + theme(plot.title=element_text(size=20, vjust=-10)) + p2[[1]] + 
        p3[[1]] + p4[[1]] + plot_layout(ncol=2, widths=c(1,8), heights=c(1,4)) & 
        theme(legend.position="none")} - {p + ggtitle("b)") + 
            theme(plot.title=element_text(size=20, vjust=-10)) + p2[[2]] + 
            p3[[2]] + p4[[2]] + plot_layout(ncol=2, widths=c(1,8), heights=c(1,4)) & 
            theme(legend.position="none")} +
            plot_layout(ncol=1)} + {legend + legend2 + plot_layout(ncol=2, width=c(1,3))} +
      plot_layout(ncol=1, heights=c(3,3,1))}
    
    if(year == "2080" & disp == "nodisp"){
      ggsave(paste0("figures/Figure_1_rev.pdf"), pall,
             width=7.5, height=9, unit="in", bg="transparent")
      #ggsave(paste0("figures/Figure_1_75_agreement.png"), pall,
      #       width=7.5, height=9, unit="in", dpi=600, bg="transparent")
    } else{
      ggsave(paste0("figures/delta_sr_", year, "_", disp, "_10CV.png"), pall,
             width=7.5, height=9, unit="in", dpi=600, bg="transparent")
      #ggsave(paste0("figures/delta_sr_", year, "_", disp, "_75agreement.png"), pall,
      #       width=7.5, height=9, unit="in", dpi=600, bg="transparent")
    }
  })
})

########################################

# Calculate relative change in probability for disp and year with taxa merged
library(dplyr); library(ggplot2); library(patchwork); library(tidyr)

# Specify year
year <- 2080

lapply(c("nodisp", "disp"), function(disp){
  # Specify taxa
  taxa <- c("Amphibian", "Ter_Bird", "Ter_Mammal")
  
  ## Create time_rcp combination
  time_rcp <- expand.grid(time=year, rcp= c("rcp26", "rcp60")) %>% 
    tidyr::unite("time_rcp", c(rcp, time))
  time_rcp <- as.vector(time_rcp$time_rcp)
  time_rcp <- c("1995", time_rcp)
  
  mean_taxa <- lapply(1:3, function(i){
    sumProb1 <- read.csv(paste0("data/", taxa[i], "_prob_GAM_", disp, ".csv.xz"))
    sumProb2 <- read.csv(paste0("data/", taxa[i], "_prob_GBM_", disp ,".csv.xz"))
    mean_taxa <- lapply(time_rcp, function(x){
      sumProb_sub1 <- sumProb1 %>% dplyr::select(c(x,y), matches(x))
      sumProb_sub2 <- sumProb2 %>% dplyr::select(c(x,y), matches(x))
      sumData <- full_join(sumProb_sub1, sumProb_sub2, by=c("x", "y")) %>% group_by(x,y)
      sumData$mean <- apply(sumData[,-c(1,2)], 1, mean, na.rm=TRUE)
      sumData <- sumData %>% dplyr::select(c(x,y,mean))
      sumData$time_rcp <- x
      return(sumData)
    })
    mean_taxa <- do.call("rbind", mean_taxa)
    mean_taxa$time_rcp <- factor(mean_taxa$time_rcp, labels=c("1995", 
                                                              paste0(year, " RCP2.6"), 
                                                              paste0(year, " RCP6.0")))
    mean_taxa$taxa <- taxa[i]
    return(mean_taxa)
  })
  mean_taxa <- bind_rows(mean_taxa)
  mean_taxa$disp <- disp
  
  mean_taxa <- tidyr::spread(mean_taxa, taxa, mean) %>% as.data.frame() %>% 
    group_by(x,y,time_rcp, disp) %>% ungroup() %>% 
    mutate(sum = rowSums(.[,5:7], na.rm=TRUE))
  
  sum_1995 <- mean_taxa %>% filter(time_rcp == 1995) %>% dplyr::select(x,y,disp,sum)
  delta_ind <- mean_taxa %>% dplyr::select(-sum) %>% left_join(sum_1995) %>%
    tidyr::gather(taxa, value, -c(x,y,time_rcp,disp, sum)) %>% 
    tidyr::spread(time_rcp, value) %>% tidyr::drop_na() %>% 
    mutate_at(c(paste0(year, " RCP2.6"), paste0(year, " RCP6.0")), funs((. - `1995`)/`sum`*100)) %>% 
    dplyr::select(-`1995`, -sum) %>% tidyr::gather(time_rcp, mean, -c(x,y,taxa,disp)) %>% 
    tidyr::spread(taxa, mean)
  delta_sum <- mean_taxa %>% dplyr::select(-c(Amphibian, Ter_Bird, Ter_Mammal)) %>% 
    tidyr::gather(taxa, value, -c(x,y,time_rcp,disp)) %>% 
    tidyr::spread(time_rcp, value) %>% tidyr::drop_na() %>% dplyr::select(-taxa) %>% 
    mutate_at(c(paste0(year, " RCP2.6"), paste0(year, " RCP6.0")), funs((. - `1995`)/`1995`*100)) %>% 
    dplyr::select(-`1995`) %>% tidyr::gather(time_rcp, sum, -c(x,y,disp))
  delta_mean <- left_join(delta_ind, delta_sum)
  
  delta_all <- delta_mean %>% dplyr::select(x,y,time_rcp, disp, sum)
  delta_x <- delta_mean %>% dplyr::select(-sum) %>% 
    tidyr::gather(taxa, mean, -c(x,y,time_rcp,disp)) %>%
    group_by(x, time_rcp, taxa, disp) %>% 
    summarise(mean=mean(mean, na.rm=T))
  delta_y <- delta_mean %>% dplyr::select(-sum) %>% 
    tidyr::gather(taxa, mean, -c(x,y,time_rcp,disp)) %>% 
    group_by(y, time_rcp, taxa, disp) %>% 
    summarise(mean=mean(mean, na.rm=T))
  
  data(outline, package="ggmap2")
  library(sf)
  outline <- sf::st_as_sf(outline)
  
  lim_map <- c(min(delta_all$sum), max(delta_all$sum))
  lim_histx <- unlist(c(delta_x %>% group_by(taxa) %>% 
                          summarise(min=min(mean, na.rm=T)) %>% 
                          ungroup() %>% summarise(sum(min, na.rm=T)), 
                        delta_x %>% group_by(taxa) %>% 
                          summarise(max=max(mean, na.rm=T)) %>% 
                          ungroup() %>% summarise(sum(max, na.rm=T))))
  lim_histy <- unlist(c(delta_y %>% group_by(taxa) %>% 
                          summarise(min=min(mean, na.rm=T)) %>% 
                          ungroup() %>% summarise(sum(min, na.rm=T)), 
                        delta_y %>% group_by(taxa) %>% 
                          summarise(max=max(mean, na.rm=T)) %>% 
                          ungroup() %>% summarise(sum(max, na.rm=T))))
  col_val <- scales::rescale(unique(c(seq(min(delta_all$sum), 0, length=5), 
                                      seq(0, max(delta_all$sum), length=5))))
  p <-  ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + theme_classic() + 
    theme(legend.position="none", axis.title = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(), 
          axis.text = element_blank())
  p2 <- lapply(c("RCP2.6", "RCP6.0"), function(rcp){
      ggplot() + geom_histogram(data=delta_x[delta_x$time_rcp == paste(year, rcp) & 
                                               delta_x$disp == disp,], 
                                aes(x=x, y=mean, fill=taxa), width=1,
                                stat="identity", position="stack", colour=NA) + 
        scale_y_continuous(limits=lim_histx, position = "right", expand=c(0,0)) + 
        scale_x_continuous(expand=c(0,0), breaks=c(15,0,-15,-30)) + 
        scale_fill_grey() + theme_classic() + 
        theme(legend.position="none", axis.title = element_blank(),
              axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), axis.text.y = element_text(size=rel(0.8)))
    })
  p3 <- lapply(c("RCP2.6", "RCP6.0"), function(rcp){
      ggplot() + geom_histogram(data=delta_y[delta_y$time_rcp == paste(year, rcp) & 
                                               delta_y$disp == disp,], 
                                aes(x=y, y=mean, fill=taxa), width=1,
                                stat="identity", position="stack", colour=NA) + 
        scale_y_reverse(limits=rev(lim_histy)) + 
        scale_fill_grey(name="", labels=c("Amphibians", "Birds", "Mammals")) + 
        coord_flip(expand=FALSE) + theme_classic() + 
        theme(legend.position="bottom", axis.title = element_blank(),
              axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
              axis.text.y = element_blank(), axis.text.x = element_text(size=rel(0.8))) + 
        guides(fill = guide_legend(direction = "vertical"))})
  #outline <- sf::st_transform(outline, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  p4 <- lapply(c("RCP2.6", "RCP6.0"), function(rcp){
      # Transform data to Robinson projection
      data <- delta_all[delta_all$time_rcp == paste(year, rcp) & 
                          delta_all$disp == disp,]
      #data <- data[,c("x","y", "sum")]
      #data <- raster::rasterFromXYZ(data)
      #raster::projection(data) <- "+proj=longlat + datum=WGS84"
      #data <- raster::projectRaster(data, crs="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
      #data <- data.frame(raster::rasterToPoints(data))
      # Make map
      ggplot() +
        geom_tile(data=data, aes(x=x, y=y, fill=sum)) + 
        geom_sf(data=outline, fill="transparent", colour="black") + 
        scale_fill_gradientn(name="", colours=rev(colorRampPalette(
          c("#00007F", "blue", "#007FFF", "cyan", 
            "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
          na.value="transparent", values=col_val, limits=lim_map) + 
        coord_sf(expand=F, xlim=c(min(data$x), max(data$x)), ylim=c(min(data$y),max(data$y)), 
                 ndiscr=0) + theme_classic() + 
        theme(axis.title = element_blank(), axis.line = element_blank(),
              axis.ticks = element_blank(), axis.text = element_blank(),
              plot.background = element_rect(fill = "transparent"), 
              legend.background = element_rect(fill = "transparent"), 
              legend.key.width=unit(2.5, "cm"), legend.position="bottom",
              legend.box.background = element_rect(fill = "transparent", colour=NA))
    })
  pall <- {{p + ggtitle("a)") + theme(plot.title=element_text(size=14, vjust=-10)) + p2[[1]] + 
      p3[[1]] + p4[[1]] + plot_layout(ncol=2, widths = c(1,8), heights=c(1,4)) & 
      theme(legend.position="none")} - {p + ggtitle("b)") + theme(plot.title=element_text(size=14, vjust=-10)) + p2[[2]] + 
          p3[[2]] + p4[[2]] + plot_layout(ncol=2, widths=c(1,8), heights=c(1,4))} +
      plot_layout(ncol=1)}
  ggsave(paste0("figures/rel_change_sr_", year, "_", disp, ".png"), pall,
         width=7.5, height=9, unit="in", dpi=600, bg="transparent")
})

########################################

## Plot coefficient of variation

# Set dispersal
library(dplyr); library(ggplot2); library(magrittr)
lapply(c("nodisp", "disp"), function(disp){
  lapply(c("2050", "2080"), function(year){
    ## Create time_rcp combination
    time_rcp <- expand.grid(time=year, rcp= c("rcp26", "rcp60")) %>% 
      tidyr::unite("time_rcp", c(rcp, time))
    time_rcp <- as.vector(time_rcp$time_rcp)
    
    # Load summed Probability 
    sumProb <- lapply(1:3, function(i){
      taxa <- c("Amphibian", "Ter_Bird", "Ter_Mammal")
      sumProb1 <- read.csv(paste0("data/", taxa[i], "_prob_GAM_", disp, ".csv.xz"))
      sumProb2 <- read.csv(paste0("data/", taxa[i], "_prob_GBM_", disp ,".csv.xz"))
      sumProb <- left_join(sumProb1, sumProb2, by=c("x","y"))
      sumProb$taxa <- c("Amphibians", "Birds", "Mammals")[i]
      return(sumProb)
    })
    sumProb %<>% bind_rows 
    
    # Summarise by taxa
    cv_taxa <- lapply(time_rcp, function(x){
      sumProb_sub <- sumProb %>% dplyr::select(c(x,y,taxa), matches(x))
      sumProb_sub$CV <- sumProb_sub %>% select(-c(x,y,taxa)) %>% apply(1, raster::cv, na.rm=TRUE)
      sumProb_sub %<>% dplyr::select(c(x,y,taxa,CV))
      sumProb_sub$time_rcp <- x
      return(sumProb_sub)
    })
    cv_taxa %<>% bind_rows()
    cv_taxa$time_rcp <- factor(cv_taxa$time_rcp, labels=c(paste0(year, " RCP2.6"), 
                                                          paste0(year, " RCP6.0")))
    cv_taxa$CV2 <- cut(cv_taxa$CV, breaks=c(0,10,20,30,40,50,60,70,80,90,100, round(max(cv_taxa$CV))))
    
    data(outline, package="ggmap2")
    p <- ggplot() +
      geom_tile(data=na.omit(cv_taxa), aes(x=x, y=y, fill=CV2)) + 
      facet_grid(taxa ~ time_rcp, scales="free") + 
      geom_polygon(data=outline, aes(x=long,y=lat, group=group), 
                   fill="transparent", colour="black") + 
      scale_fill_manual(name="%", values=colorRampPalette(
        c("white", "#00007F", "blue", "#007FFF", "cyan", 
          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(11),
        labels=c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", ">100"),
        na.value="transparent") + 
      theme_classic() + theme(axis.title = element_blank(),axis.line = element_blank(),
                              axis.ticks = element_blank(), axis.text = element_blank(),
                              panel.grid = element_blank(), 
                              strip.background= element_blank(),
                              panel.background = element_rect(fill = "transparent"),
                              plot.background = element_rect(fill = "transparent"),
                              strip.text = element_text(size=12, face="bold", colour="black"),
                              legend.background = element_rect(fill = "transparent"),
                              legend.box.background = element_rect(fill = "transparent", colour=NA)) + 
      coord_quickmap(xlim=c(-180,180), ylim=c(-60,85), expand=FALSE)
    ggsave(file=paste0("figures/ensemble_sumProb_CV_", year, "_", disp, "_rev.png"), p,
           width=9, height=6, unit="in", dpi=300, bg="transparent")
  })
})

########################################

## Plot different dispersal scenarios

# Set dispersal
library(dplyr); library(ggplot2); library(patchwork)
lapply(c("2050", "2080"), function(year){
  mean_taxa <- lapply(1:3, function(i){
    taxa <- c("Amphibian", "Ter_Bird", "Ter_Mammal")
    mean <- lapply(c("disp1", "disp2", "disp3", "disp4"), function(disp){
      sumProb1 <- read.csv(paste0("data/", taxa[i], "_prob_GAM_", disp, ".csv.xz"))
      sumProb2 <- read.csv(paste0("data/", taxa[i], "_prob_GBM_", disp ,".csv.xz"))
      ## Create time_rcp combination
      time_rcp <- expand.grid(time=year, rcp= c("rcp26", "rcp60")) %>% 
        tidyr::unite("time_rcp", c(rcp, time))
      time_rcp <- as.vector(time_rcp$time_rcp)
      mean_taxa <- lapply(time_rcp, function(x){
        sumProb_sub1 <- sumProb1 %>% select(c(x,y), matches(x))
        sumProb_sub2 <- sumProb2 %>% select(c(x,y), matches(x))
        sumData <- full_join(sumProb_sub1, sumProb_sub2, by=c("x", "y")) %>% group_by(x,y)
        sumData$mean <- apply(sumData[,-c(1,2)], 1, mean, na.rm=TRUE)
        sumData <- sumData %>% select(c(x,y,mean))
        sumData$time_rcp <- x
        return(sumData)
      })
      mean_taxa <- do.call("rbind", mean_taxa)
      mean_taxa$time_rcp <- factor(mean_taxa$time_rcp, labels=c(paste0(year, " RCP2.6"), 
                                                                paste0(year, " RCP6.0")))
      mean_taxa$disp <- disp
      return(mean_taxa)
    })
    mean <- do.call("rbind", mean)
    mean$taxa <- taxa[i]
    return(mean)
  })  
  mean_taxa <- do.call("rbind", mean_taxa)
  
  # Summarise Taxa
  mean <- mean_taxa %>% group_by(x,y,time_rcp,disp) %>%
    summarise(mean = sum(mean))
  
  # Plot map
  data(outline, package="ggmap2")
  mean$disp <- as.factor(mean$disp)
  mean$disp <- factor(mean$disp,levels(mean$disp)[c(5,1:4)], 
                      labels=c("d/4", "d/2", "d", "2*d"))
  p_mean <- ggplot() +
    geom_tile(data=mean, aes(x=x, y=y, fill=mean)) + 
    facet_grid(disp ~ time_rcp) + 
    geom_polygon(data=outline, aes(x=long,y=lat, group=group), 
                 fill="transparent", colour="black") + 
    scale_fill_gradientn(name="SR", limits=c(0, max(mean$mean)),
                         colours=colorRampPalette(
                           c("white", "#00007F", "blue", "#007FFF", "cyan", 
                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(255), 
                         na.value="transparent", values=scales::rescale(seq(0,1, length.out=10))) + 
    theme_classic() + theme(axis.title = element_blank(),axis.line = element_blank(),
                            axis.ticks = element_blank(), axis.text = element_blank(),
                            panel.grid = element_blank(), 
                            strip.background= element_blank(),
                            panel.background = element_rect(fill = "transparent"),
                            plot.background = element_rect(fill = "transparent"),
                            legend.background = element_rect(fill = "transparent"),
                            legend.box.background = element_rect(fill = "transparent", colour=NA),
                            strip.text.x = element_text(size=12, face="bold", colour="black"),
                            strip.text.y = element_text(size=10, face="bold", colour="black")) + 
    coord_quickmap(xlim=c(-180,180), ylim=c(-60,85), expand=FALSE)
  
  # Calculate CV among different dispersal scenarios
  cv <- mean %>% filter(disp != "nodisp") %>% select(-disp) %>% 
    group_by(x, y, time_rcp) %>% 
    summarise(cv=raster::cv(mean, na.rm=TRUE))
  cv$cv <- cut(cv$cv, breaks=c(0,20,40,60,80,100, round(max(cv$cv))))
  
  p_cv <- ggplot() +
    geom_tile(data=cv, aes(x=x, y=y, fill=cv)) + 
    facet_wrap(~ time_rcp, nrow=1) + 
    geom_polygon(data=outline, aes(x=long,y=lat, group=group), 
                 fill="transparent", colour="black") + 
    scale_fill_manual(name="%", values=c("white", "#00007F", "#007FFF",
        "#7FFF7F", "#FF7F00", "#7F0000"), labels=c("0-20", "20-40", "40-60", "60-80", "80-100", ">100")) + 
    theme_classic() + theme(axis.title = element_blank(),axis.line = element_blank(),
                            axis.ticks = element_blank(), axis.text = element_blank(),
                            panel.grid = element_blank(), 
                            strip.background= element_blank(),
                            panel.background = element_rect(fill = "transparent"),
                            plot.background = element_rect(fill = "transparent"),
                            legend.background = element_rect(fill = "transparent"),
                            legend.box.background = element_rect(fill = "transparent", colour=NA),
                            strip.text = element_blank()) + 
    coord_quickmap(xlim=c(-180,180), ylim=c(-60,85), expand=FALSE)
  pall <- p_mean + p_cv + plot_layout(ncol=1, width = c(1,1), height=c(4,1))
  ggsave(file=paste0("figures/sumProb_disp_", year, "_rev.png"), pall,
         width=8, height=8, unit="in", dpi=600, bg="transparent")
})
