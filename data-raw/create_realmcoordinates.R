#' ## Create Land coordinates with terrestrial realm information

#' Land coordinates from ISIMIP
library(raster); library(ggplot2)
land <- raster("extdata/ISIMIP2b_landseamask_generic.nc4")
df_land <- data.frame(rasterToPoints(land))
colnames(df_land) <- c("x", "y", "LSM")
ggplot()+geom_raster(data=df_land,aes(x=x,y=y, fill=LSM))

#' Read Terrestrial ecoregions of the world shapefile
teow <- rgdal::readOGR("extdata/wwf_terr_ecos.shp")
teow_realm <- gUnaryUnion(teow, teow$REALM)
#' which was downloaded from: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
plot(teow_realm)

#' Alke suggested to use: CMEC Zoogeographic Realms and Regions
#' which can be downloaded here: http://macroecology.ku.dk/resources/wallace#Gis
zoo_realm <- rgdal::readOGR("extdata/newRealms.shp")
plot(zoo_realm, col=zoo_realm$Realm)

#' Rasterize CMEC according to ISIMIP2b landseamask
#' using all cells
line <- as(zoo_realm, "SpatialLines")
r_line <- raster::rasterize(line, land, background=NA, na.rm=TRUE)
r_poly <- raster::rasterize(zoo_realm, land, background=NA, na.rm=TRUE)
r_realm <- raster::merge(r_line, r_poly)

#'using the center of the cell
r_realm <- raster::rasterize(zoo_realm, land, background=NA, na.rm=TRUE)

#' Turn raster into a dataframe
realm_coordinates<- as.data.frame(rasterToPoints(r_realm))
realm_coordinates <- realm_coordinates[,c("x", "y", "layer")]
colnames(realm_coordinates)[3] <- "Realm"
realm_coordinates$Realm <- as.factor(realm_coordinates$Realm)

# Merge land and realm
library(dplyr)
realm_coordinates <- left_join(df_land, realm_coordinates, by=c("x","y"))

#' Plot realm coordinates
ggplot()+geom_raster(data=df_final, aes(x=x,y=y, fill=Realm))

# Add area to realm_coordinates
data(landseamask_generic, package="rISIMIP")
area <- as.data.frame(raster::rasterToPoints(raster::area(landseamask_generic, na.rm=TRUE)))
colnames(area) <- c("x", "y", "area")
realm_coordinates <- dplyr::full_join(realm_coordinates, area, by=c("x","y"))

#' Save dataframe to file
write.csv2(realm_coordinates, "data/realm_coordinates.csv", row.names=FALSE)

