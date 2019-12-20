# Top land-use rev

# Specify path of file directory
filedir <- "/media/mbiber/BioScen1point5_Data"

# Choose one of totals, 5crops, 15crops
lu_type <- "totals"

# Get summarised landuse data
crops <- list.files(paste0(filedir, "/ISIMIP2b/DerivedInputData/landuse/global/"), 
                    pattern=glob2rx(paste0("*", lu_type, "*.csv$")), full.names=T)

# Read files into dataframe
crops <- lapply(crops, function(x){
  data <- read.csv(x)
  data$year <- strsplit(strsplit(basename(x), split="_")[[1]][4], split="[.]")[[1]][1]
  data$scenario <- strsplit(basename(x), split="_")[[1]][1]
  data$model <- strsplit(basename(x), split="_")[[1]][2]
  return(data)
})
crops <- do.call(plyr::rbind.fill, crops)

#Calculate area of each cell in km2
data(landseamask_generic, package="rISIMIP")
isimip_area <- raster::area(landseamask_generic, na.rm=TRUE) # km2
isimip_area <- as.data.frame(raster::rasterToPoints(isimip_area))
colnames(isimip_area) <- c("x", "y", "area")
sum(isimip_area$area)

# Add area to crops dataframe, remove 2005soc data
library(dplyr)
crops %<>% left_join(isimip_area) %>% filter(scenario != "2005soc")

landuse_all <- crops %>% dplyr::select(-c(area, matches("total"))) 
readr::write_csv(landuse_all, "data/landuse_all.csv.xz")
