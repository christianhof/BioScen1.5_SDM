#' Code to summarise the model output

# Load dplyr package
rm(list=ls()); gc()
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(snowfall)

# Set file directory
filedir <- "I:/"

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird", "Reptile")
i <- 3

# Model type
k <- 4; model_type <- c("GAM", "GBM", "MaxEnt", "RF")[k]

# Time steps
if(taxa[i]=="Reptile"){timesteps <- c(1995, 2050, 2080)} else{
  timesteps <- c(1845, 1990, 1995, 2009, 2010, 2020, 2026, 2032, 2048, 2050, 
                 2052, 2056, 2080, 2100, 2150, 2200, 2250)
}

# File directory of results
predictsPaths <- sapply(timesteps, function(x){
  paste0(filedir, "/", taxa[i],  "_", model_type, "_predictions_", x)
})

# Output directory
sum_predPath <- paste0(filedir, "/", taxa[i], "_", model_type, "_predictions")
if(!dir.exists(sum_predPath)){dir.create(sum_predPath)}

# List all files
Modelfiles <- lapply(predictsPaths, function(x) list.files(x, full.names=TRUE))
Modelfiles <- do.call("c", Modelfiles)

# Check one output
Model1File <- read.csv(list.files(predictsPaths[[1]], full.names=TRUE)[1])
#head(Model1File)


# Aggregate the different AUC values from the 10 iterations per species
# and filter by AUC > 0.7
if(taxa[i]=="Reptile"){
  # Extract species names 
  AUC_data <- lapply(c("GAM", "GBM"), function(model_type){
    read.csv(paste0(filedir, "/AUCvalues_All_", 
                    model_type, "_", taxa[i], ".csv.xz"))})
  AUC_data <- do.call(rbind, AUC_data)
  
  AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
    dplyr::summarise(mean = mean(AUC, na.rm=T)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
    group_by(Species, taxa) %>% dplyr::summarise(n = n()) %>% filter(n == 2)
  } else{
  # Extract species names 
  AUC_data <- lapply(c("GAM", "GBM", "MaxEnt", "RF"), function(model_type){
    read.csv(paste0(filedir, "/AUCvalues_All_", 
                    model_type, "_", taxa[i], ".csv.xz"))})
  AUC_data <- do.call(rbind, AUC_data)
  AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
    dplyr::summarise(mean = mean(AUC, na.rm=T)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
    group_by(Species, taxa) %>% dplyr::summarise(n = n()) %>% filter(n == 4)
}
spNames <- unique(AUC_sum$Species)
length(spNames)

# Get missing names
names_mis <- lapply(spNames, function(x){
  if(!file.exists(paste0(sum_predPath, "/", x, "_", model_type, "_predict.csv.xz"))){
    return(x)
  }
})
names_mis <- unlist(Filter(Negate(is.null), names_mis))
# Remove Ursus maritimus
names_mis <- names_mis[!names_mis %in% "Ursus_maritimus"]
length(names_mis)

#' Check of model prediction files

# List all prediction files
#files <- unlist(lapply(dir(filedir, pattern=paste0(model_type, "_predictions_"), full.names = T), 
#                   function(x){list.files(x, full.names=T)}))
#length(files)

##Get all files for one species
#files <- unlist(lapply(names_mis, function(species){files[grepl(files, pattern=species)]}))
#length(files)
#head(files)

# Check for corrupt files
#snowfall::sfInit(parallel=TRUE, cpus=ceiling(0.25*parallel::detectCores()))
#corrupt_files <- snowfall::sfLapply(files, function(x){
#  data <- tryCatch(readr::read_csv(x), error=function(e) e) #The prediction files
#  if(inherits(data, "error")){
#    return(x)
#  }
#}); snowfall::sfStop()

#corrupt_files <- unlist(Filter(Negate(is.null), corrupt_files))
#length(corrupt_files)
#file.remove(corrupt_files) # Remove corrupt files

#' Loop through all species names and save summarized prediction
n <- 10
sfInit(parallel=TRUE, cpus=n)
sfLibrary(dplyr); sfLibrary(readr); sfLibrary(tidyr)
sfExport(list=c("Modelfiles", "sum_predPath", "model_type", "timesteps")) 

#Run in parallel
lapply(names_mis, function(species){
  ##Get all files for one species
  spFiles <- Modelfiles[grepl(Modelfiles, pattern=paste0(species, "_"))]
  
  avg_data <- lapply(timesteps, function(z){
    print(paste(species, z))
    ##Import data for all species separated by timeframe
    Files <- grep(spFiles, pattern=z, value=TRUE)
    
    # Read data and merge into one dataframe
    data <- readr::read_csv(Files)
    colnames(data) <- gsub("[.]", "-", colnames(data))
    #data <- lapply(Files, function(x){read.csv(paste0(x))})
    #data <- do.call("bind_rows", data)
    
    ## Create model_rcp combination
    model <- c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5")
    if(z == 1845){
      rcp <- "piControl"
      model_rcp <- expand.grid(model=model, rcp=rcp)
      model_rcp <- tidyr::unite(model_rcp, "model_rcp", c(model,rcp))
      model_rcp <- as.vector(model_rcp$model_rcp)
    } else if(z == 1990){
      rcp <- "historical"
      model_rcp <- expand.grid(model=model, rcp=rcp)
      model_rcp <- tidyr::unite(model_rcp, "model_rcp", c(model,rcp))
      model_rcp <- as.vector(model_rcp$model_rcp)
    } else if(z == 1995){
      model_rcp <- "EWEMBI"
    } else if(model_type %in% c("GAM", "GBM") & z %in% c(2050, 2080)){
      rcp <- c("rcp26", "rcp60", "rcp85")
      model_rcp <- expand.grid(model=model, rcp=rcp)
      model_rcp <- tidyr::unite(model_rcp, "model_rcp", c(model,rcp))
      model_rcp <- as.vector(model_rcp$model_rcp)
    } else if(z < 2100){
      rcp <- c("rcp26", "rcp60")
      model_rcp <- expand.grid(model=model, rcp=rcp)
      model_rcp <- tidyr::unite(model_rcp, "model_rcp", c(model,rcp))
      model_rcp <- as.vector(model_rcp$model_rcp)
    } else{
      model_rcp <- c("HadGEM2-ES_rcp26", "IPSL-CM5A-LR_rcp26", "MIROC5_rcp26")
    }
    
    period <- z
    if(colnames(data)[1]!="x"){
      n <- ncol(data)-2
      if(period == 1995){
        y <- 1
        colnames(data)[3:length(colnames(data))] <- paste0(model_type, "_EWEMBI_1995_block_", seq(1:(n/y)))
        colnames(data)[1] <- "x"
        colnames(data)[2] <- "y"
      } else if(period == 1845){
        y <- 4
        colnames(data)[3:length(colnames(data))] <- c(paste0(model_type, "_GFDL-ESM2M_piControl_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_HadGEM2-ES_piControl_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_IPSL-CM5A-LR_piControl_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_MIROC5_piControl_block_", seq(1:(n/y))))
        colnames(data)[1] <- "x"
        colnames(data)[2] <- "y"
      } else if(period == 1990){
        y <- 4
        colnames(data)[3:length(colnames(data))] <- c(paste0(model_type, "_GFDL-ESM2M_historical_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_HadGEM2-ES_historical_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_IPSL-CM5A-LR_historical_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_MIROC5_historical_block_", seq(1:(n/y))))
        colnames(data)[1] <- "x"
        colnames(data)[2] <- "y"
      }else if(period <= 2080){
        y <- 8
        colnames(data)[3:length(colnames(data))] <- c(paste0(model_type, "_GFDL-ESM2M_rcp26_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_GFDL-ESM2M_rcp60_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_HadGEM2-ES_rcp26_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_HadGEM2-ES_rcp60_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_IPSL-CM5A-LR_rcp26_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_IPSL-CM5A-LR_rcp60_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_MIROC5_rcp26_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_MIROC5_rcp60_block_", seq(1:(n/y))))
        colnames(data)[1] <- "x"
        colnames(data)[2] <- "y"
      } else{
        y <- 3 
        colnames(data)[3:length(colnames(data))] <- c(paste0(model_type, "_HadGEM2-ES_rcp26_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_IPSL-CM5A-LR_rcp26_block_", seq(1:(n/y))),
                                                    paste0(model_type, "_MIROC5_rcp26_block_", seq(1:(n/y))))
        colnames(data)[1] <- "x"
        colnames(data)[2] <- "y"
      }
      readr::write_csv(x=data, path=Files)
    }
    
    ## Select all data for one model_rcp combination
    ## And calculate average among blocks and PAs
    avg_data <- lapply(model_rcp, FUN=function(w){
      sub_data <- data %>% dplyr::select(x, y, contains(w)) %>% drop_na() %>%
        tidyr::gather(var, value, -c(x,y)) %>% 
        dplyr::group_by(x,y) %>% dplyr::summarise(mean=round(mean(value, na.rm=TRUE), 3))
      colnames(sub_data) <- c("x", "y", paste0(w, "_", z))
      return(sub_data)
    })
    avg_data <- Reduce(function(x, y) full_join(x, y, by=c("x","y")), avg_data)
    if(nrow(avg_data)==0){
      avg_data[1,] <- NA
      avg_data$x <- as.numeric(avg_data$x)
      avg_data$y <- as.numeric(avg_data$y)
    }
    gc()
    return(avg_data)
  })
  
  # Combine data
  avg_data <- Reduce(function(x,y) full_join(x, y, by=c("x", "y")), avg_data)
  
  # Remove entries where all data is 0
  avg_data <- avg_data[!!rowSums(abs(avg_data[-c(1,2)]), na.rm=T),]
  
  # Remove entries where all data are NA
  avg_data <- avg_data[rowSums(is.na(avg_data)) != ncol(avg_data), ]
  
  # Save data to csv file
  readr::write_csv(avg_data, path=paste0(sum_predPath, "/",  species, "_", 
                                         model_type, "_predict.csv.xz"))
  return(NULL)
})
sfStop() # Close the cluster
#system('shutdown -s')

# Test output for one species
test <- read.csv(list.files(sum_predPath, full.names=T)[[1]])
head(test)

# Plot prediction
#ggplot(data=test, aes(x=x, y=y, fill=MIROC5_rcp26_2020)) + geom_raster()

########################################

#' Code to save individual model predictions in one NetCDF file

# Load packages
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ncdf4)

# Set file directory
#filedir <- "/scratch/home/mbiber/data" # shinichi
#filedir <- "/bigdata_local/mbiber" # ceremony - Mammals
#filedir <- "/home/mbiber/data" # ceremony - Birds
filedir <- "/bigdata/mbiber/data" # Darkstar
#filedir <- "/scratch/mbiber/data"
#filedir <- "H:/"

# Set taxa
taxa <- c("Amphibian", "Mammal", "Bird")
i <- 3

# Output directory
sum_predPath <- paste0(filedir, "/OutputData/biodiversity/")
if(!dir.exists(sum_predPath)){dir.create(sum_predPath, recursive=T)}

# Read data
AUC_data <- lapply(c("GAM", "GBM", "MaxEnt", "RF"), function(model_type){
  read.csv(paste0(filedir, "/AUCvalues_All_", 
                  model_type, "_", taxa[i], ".csv.xz"))})
AUC_data <- do.call(rbind, AUC_data)

# Aggregate the different AUC values from the 10 iterations per species
# and filter by AUC > 0.7
AUC_sum <- AUC_data %>% group_by(Species, taxa, model_type) %>% 
  summarise(mean = mean(AUC, na.rm=T)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
  group_by(Species, taxa) %>% summarise(n = n()) %>% filter(n == 4)

spNames <- unique(AUC_sum$Species); rm(AUC_data, AUC_sum)
# Remove Ursus maritimus
spNames <- spNames[!spNames == "Ursus_maritimus"]
length(spNames)

# Define model/rcp combinations
model_type <- "RF" #c("GAM", "GBM", "MaxEnt", "RF")
df <- expand.grid(model_type=model_type, 
                  rcp=c("historical", "piControl", "rcp26", "rcp60"))
df <- rbind(expand.grid(model_type=model_type, rcp="1995"), df)
df <- rbind(df, expand.grid(model_type=model_type, rcp="2100rcp26"))
df

#' Loop through all species names and save NetCDF file

# Run code for each model output file
rcp_mod <- 4

#Define years
if(df$rcp[rcp_mod] == "piControl"){
  year <- 1845
} else if(df$rcp[rcp_mod] == "1995"){
  year <- 1995
} else if(df$rcp[rcp_mod] == "historical"){
  year <- 1990
} else if(df$rcp[rcp_mod] %in% c("rcp26", "rcp60")){
  year <- c(2009, 2010, 2020, 2026, 2032, 2048, 2050, 2052, 2056, 2080)
} else if(df$rcp[rcp_mod] == "2100rcp26"){
  year <- c(2100, 2150, 2200, 2250)
  df$rcp[rcp_mod] <- "rcp26"
}

# Define predicts paths according to years
predictsPaths <- sapply(year, function(x){paste0(filedir, "/", taxa[i], "_", 
                                                 df$model_type[rcp_mod], "_predictions_", x)})
if(df$rcp[rcp_mod] == "1995"){
  model="EWEMBI"
  filename <- paste0(filedir, "/OutputData/biodiversity/bioscen1.5-sdm-", 
                     tolower(df$model_type[rcp_mod]),
                     "_ewembi_nobc_hist_nosoc_co2_", tolower(taxa[i]), 
                     "prob_global_30year-mean_1995_1995.nc4")
} else{
  if(any(year >= 2100)){
    model=c("MIROC5", "HadGEM2.ES", "IPSL.CM5A.LR")
    filename <- sapply(model, function(m){paste0(filedir, "/OutputData/biodiversity/bioscen1.5-sdm-", tolower(df$model_type[rcp_mod]), "_",
                                                 gsub("[.]", "-", tolower(m)), "_ewembi_2100", 
                                                 tolower(df$rcp[rcp_mod]), "_nosoc_co2_", 
                                                 tolower(taxa[i]), "prob_global_30year-mean_", 
                                                 min(year), "_", max(year),".nc4")})
  } else{
    model=c("MIROC5", "HadGEM2.ES", "IPSL.CM5A.LR", "GFDL.ESM2M")
    filename <- sapply(model, function(m){paste0(filedir, "/OutputData/biodiversity/bioscen1.5-sdm-", tolower(df$model_type[rcp_mod]), "_",
                                                 gsub("[.]", "-", tolower(m)), "_ewembi_", 
                                                 tolower(df$rcp[rcp_mod]), "_nosoc_co2_", 
                                                 tolower(taxa[i]), "prob_global_30year-mean_", 
                                                 min(year), "_", max(year),".nc4")})
  }
}
length(which(!file.exists(filename)))

if(length(which(!file.exists(filename))) > 0){
  #Define the dimensions
  dimX = ncdim_def(name="lon", units="degrees", vals=seq(-179.75, 179.75, length = 720))
  dimY = ncdim_def(name="lat", units="degrees", vals=seq(89.75, -89.75, length = 360))
  dimT = ncdim_def(name="time", units="years since 1661-1-1 00:00:00", 
                   vals=c(year-1661), calendar="proleptic_gregorian")
  
  # Define data for NetCDF file
  vard <- lapply(spNames, function(name){
    ncvar_def(as.character(name), "Probability of occurrence per cell", 
              list(dimX,dimY,dimT), 1.e+20, prec="double", compression=9)})
  
  # Create the NetCDF files
  filename <- filename[which(!file.exists(filename))]
  model <- model[which(!file.exists(filename))]
  lapply(filename, function(x){
    nc <- nc_create(x, vard)
    ncatt_put(nc, varid=0, attname="contact", attval="Matthias Biber <matthias.biber@tum.de>")
    ncatt_put(nc, varid=0, attname="institution", attval="Technical University Munich (Germany)")
    nc_close(nc)
  })
  
  # Individually write data for every species
  for(j in 6770:length(spNames)){
    library(readr)
    print(j)
    # Get data
    files <- lapply(predictsPaths,function(x){list.files(x, pattern=paste0(as.character(spNames[j]), "_"), 
                                                         full.names=T)})
    if(length(files) == 1){
      data <- readr::read_csv(files[[1]])
      data$PA <- 1
      data$year <- year
    } else{
      data <- lapply(1:length(files), function(x){
        data <- readr::read_csv(files[[x]])
        data$PA <- 1
        data$year <- year[x]
        return(data)
      })
      data <- do.call(plyr::rbind.fill, data)
    }
    
    ## Select all data for one model_rcp combination
    ## And calculate average among blocks and PAs
    # Spread years
    for(k in 1:length(filename)){
      nc <- nc_open(filename[k], write=T)
      
      data_sub <- data %>% group_by(x, y, PA, year) %>% 
        select(x, y, PA, year,matches(paste(model[k],df$rcp[rcp_mod], sep="_"))) %>% 
        tidyr::gather(var, value, -c(x,y,year,PA)) %>% group_by(x,y,year) %>% 
        summarise(mean=mean(value, na.rm=TRUE)) %>% tidyr::spread(year, mean)
      
      #Expand dataframe with NAs
      df_spat <- expand.grid(x=seq(-179.75, 179.75, length = 720), 
                             y=seq(89.75, -89.75, length = 360))
      data_sub <- left_join(df_spat, data_sub) %>% select(-x, -y); rm(df_spat)
      
      # Turn data into array
      data_sub <- array(unlist(data_sub),dim=c(720, 360, ncol(data_sub)), 
                        dimnames=list(NULL, NULL, names(data_sub)))
      
      # Write data to the NetCDF file
      ncvar_put(nc, vard[[j]], data_sub, start=c(1,1,1), count=c(-1,-1,-1))
      
      # Close your new file to finish writing
      nc_close(nc)
    }
  }
}
sfStop()
#system('shutdown -s')

########################################

# Test NetCDF file
library(raster); library(ncdf4)
filedir <- "C:/Users/admin/Documents/"
i <- list.files(filedir, pattern = ".nc4", full.names=T)
(files <- list.files(paste0(filedir, "OutputData/biodiversity/"), pattern="gam.*historical", full.names=T))
#i <- list.files(paste0(filedir, "OutputData/biodiversity/"), pattern="maxent_IPSL-CM5A-LR_rcp60", full.names=T)
#i <- files[1]
#j <- "Yunganastes_mercedesae"
for(i in files){
  par(mfrow=c(2,2))
  for(j in c("Acanthixalus_spinosus", "Limnodynastes_convexiusculus", 
             "Xenorhina_rostrata", "Yunganastes_mercedesae")){
    nc <- nc_open(i)
    v <- nc$var[[nc$nvars]]
    #nc_close(nc)
    #test <- stack(i, varname=v$name)
    test <- stack(i, varname=j)
    print(nlayers(test))
    plot(test[[sample(1:nlayers(test),1)]])
  }
}