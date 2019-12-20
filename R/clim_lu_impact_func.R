clim_lu_impact <- function(x){
  print(x)
  
  # Read GAM and GBM data
  gam <- read.csv(grep(files_GAM, pattern=paste0(x, "_"), value=T))
  gbm <- read.csv(grep(files_GBM, pattern=paste0(x, "_"), value=T))
  gam$model <- "GAM"
  gbm$model <- "GBM"
  spdata <- bind_rows(gam, gbm)
  
  if(length(spdata) != 0){
    # Turn data into long format
    library(dplyr)
    library(magrittr)
    spdata <- select(spdata, x,y,areaKM2, presence, model,
                     EWEMBI_1995:MIROC5_rcp60_2080) %>% 
      tidyr::gather(var, value, -c(x,y,areaKM2, presence, model)) 
    spdata <- spdata %>% tidyr::separate(col=var, c("GCM", "scenario", "year"), 
                                         sep="_", remove=T, extra="drop", fill="left")
    spdata$GCM[spdata$scenario == "EWEMBI"] <- "EWEMBI"
    spdata$scenario[spdata$GCM == "EWEMBI"] <- NA
    spdata <- tidyr::unite(spdata, scenario, year, col="time_rcp", sep=" ")
    spdata$time_rcp[spdata$time_rcp == "NA 1995"] <- "1995"
    
    # Get separate data.frame for EWEMBI Data
    ewembi <- spdata %>% filter(GCM == "EWEMBI") %>% select(-c(GCM, time_rcp))
    ewembi$`1995` <- ewembi$value
    ewembi <- ewembi %>% select(-value)
    head(ewembi)
    
    # Change format back to wide
    sumData <- spdata %>% filter(GCM != "EWEMBI") %>% 
      tidyr::spread(time_rcp, value)
    
    # Add EWEMBI to sumData
    sumData <- left_join(sumData, ewembi)
    
    #Calculate threshold
    #sumData$thres <- weighted.mean(sumData$`1995`, sumData$areaKM2, na.rm=T)
    
    # Mutate cells and identify cells where climate suitability is decreasing (TRUE)
    sum_disp_thres <- sumData %>% 
      mutate_at(vars(`rcp26 2080`:`rcp60 2080`), 
                funs(. < weighted.mean(`1995`, `areaKM2`, na.rm=T)))
    sum_disp_thres$disp <- 1
    sum_disp_thres$thres <- 1
    
    sum_nodisp_thres <- sumData %>% filter(presence == 1) %>% 
      mutate_at(vars(`rcp26 2080`:`rcp60 2080`), 
                funs(. < weighted.mean(`1995`, `areaKM2`, na.rm=T)))
    sum_nodisp_thres$disp <- 0
    sum_nodisp_thres$thres <- 1

    sum_disp_nothres <- sumData %>% 
      mutate_at(vars(`rcp26 2080`:`rcp60 2080`), funs(. < `1995`))
    sum_disp_nothres$disp <- 1
    sum_disp_nothres$thres <- 0
    
    sum_nodisp_nothres <- sumData %>% filter(presence == 1) %>% 
      mutate_at(vars(`rcp26 2080`:`rcp60 2080`), funs(. < `1995`))
    sum_nodisp_nothres$disp <- 0
    sum_nodisp_nothres$thres <- 0
    sum_disp_all <- bind_rows(sum_nodisp_nothres, sum_disp_nothres, 
                              sum_nodisp_thres, sum_disp_thres)
    
    # Turn species data into long format
    sum_disp_all %<>% select(-c(`1995`, presence)) %>% 
      tidyr::gather(time_rcp, value, -c(x,y,areaKM2, model, GCM, thres, disp))

    # Split time_rcp into scenario and year
    sum_disp_all %<>% tidyr::separate(col=time_rcp, c("scenario", "year"), sep=" ")
    
    # Merge species data with landuse data
    sum_disp_all %<>% left_join(landuse_change, 
                                by=c("x","y","GCM", "model", "scenario", "year"))
    sum_disp_all$biofuel_cropland[is.na(sum_disp_all$biofuel_cropland)] <- FALSE
    sum_disp_all$cropland[is.na(sum_disp_all$cropland)] <- FALSE
    sum_disp_all$pastures[is.na(sum_disp_all$pastures)] <- FALSE
    
    #ggplot() + geom_raster(sum_disp_all, aes(x,y,fill=biofuel_cropland))
    
    # Calculate area that is affected by combination of climate and LU
    sum_disp_all$threat[sum_disp_all$value == 1] <- 1
    sum_disp_all$threat[sum_disp_all$biofuel_cropland == 1] <- 2
    sum_disp_all$threat[sum_disp_all$cropland == 1] <- 3
    sum_disp_all$threat[sum_disp_all$pastures == 1] <- 4
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$biofuel_cropland == 1] <- 5
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$cropland == 1] <- 6
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$pastures == 1] <- 7
    sum_disp_all$threat[sum_disp_all$biofuel_cropland == 1 & sum_disp_all$cropland == 1] <- 8
    sum_disp_all$threat[sum_disp_all$biofuel_cropland == 1 & sum_disp_all$pastures == 1] <- 9
    sum_disp_all$threat[sum_disp_all$cropland == 1 & sum_disp_all$pastures == 1] <- 10
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$biofuel_cropland == 1 & 
                                sum_disp_all$cropland == 1] <- 11
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$biofuel_cropland == 1 & 
                                sum_disp_all$pastures == 1] <- 12
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$cropland == 1 & 
                                sum_disp_all$pastures == 1] <- 13
    sum_disp_all$threat[sum_disp_all$biofuel_cropland == 1 & sum_disp_all$cropland == 1 & 
                                sum_disp_all$pastures == 1] <- 14
    sum_disp_all$threat[sum_disp_all$value == 1 & sum_disp_all$biofuel_cropland == 1 & 
                                sum_disp_all$cropland == 1 & sum_disp_all$pastures == 1] <- 15
    
    sum_disp_all %<>% ungroup() %>% select(-c(x,y,value,biofuel_cropland, cropland, pastures)) %>%
      group_by(scenario, year, model, GCM, threat, disp, thres) %>% 
      summarise(sum=sum(areaKM2,na.rm=T))
    sum_disp_all$species <- x
    return(sum_disp_all)
  }
}