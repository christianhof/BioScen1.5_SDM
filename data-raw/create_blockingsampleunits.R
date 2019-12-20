#' # Create blocking units

library(raster)
library(rgdal)

#' ## Set up 0.5 degree Raster and reshape the ecoregion data

# Set coordinates for regional extent
xmin <- -180
xmax <- 180
ymin <- -90
ymax <- 90

# change these values to correct resolution
# Create 0.44 x 0.440002 degree raster layer grid and set 
# coordinates to same as climate raster (suffix: l=0.5/1 degree, s=2.5')
r.grid.l <- raster(nrows=360, ncols=720, xmn=xmin, xmx=xmax, ymn=ymin, 
                   ymx=ymax,crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
vals <- 1:ncell(r.grid.l) #create vector of numbers
r.grid.l <- setValues(r.grid.l, vals) #fill grid squares with numerical value to create label
r.grid.s <- disaggregate(r.grid.l, fact=c(10,10),method="") #disaggregate so that each 2.5' cell has number associating it to a 1 degree cell

# Create grid for dividing large ecoregions
r.grid.20 <- raster(nrows=360/10, ncols=720/10, xmn=xmin, xmx=xmax, ymn=ymin, 
                    ymx=ymax,crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
vals <- 1:ncell(r.grid.20) #create vector of numbers
r.grid.20.l <- setValues(r.grid.20, vals) #fill grid squares with numerical value to create label
r.grid.20.s <- disaggregate(r.grid.20.l, fact=c(10,10),method="") #disaggregate so that each 2.5' cell has number associating it to a 1 degree cell
r.grid.20.s <-crop(r.grid.20.s,c(xmin,xmax,ymin,ymax))

# Import ecoregion data
i <- "extdata/wwf_terr_ecos.shp"
ecoReg <- readOGR(i, "wwf_terr_ecos")
ecoReg.num <- ecoReg[8] #ecoregion id
ecoReg.area <- ecoReg[2] #ecoregion area

# Rasterize and convert into sampling units with max size of 10*res

# Can rasterize also be called on multiple things, to output stack???
# Would save re-running things 3 times.
raster.feature <- rasterize(ecoReg.num, r.grid.s)
raster.feature.l <- aggregate(raster.feature, 10, modal, progress='text') #if it overlaps with ecor 1, make it ecoregion 1
raster.ecoReg.num <- rasterize(ecoReg.num, r.grid.s, field=names(ecoReg.num))
raster.ecoReg.num.l <- aggregate(raster.ecoReg.num, 10, modal, progress='text')
raster.ecoReg.area <- rasterize(ecoReg.area, r.grid.s, field=names(ecoReg.area))
raster.ecoReg.area.l <- aggregate(raster.ecoReg.area, 10, modal, progress='text')
r.feat.reg <- stack(raster.feature.l,raster.ecoReg.num.l,r.grid.20.s) #each 0.5 accociated with ecoregion and larger grid cell
reg.id <- function(x,na.rm){as.numeric(paste(x[2],x[3],sep="."))}
sample.unit.id <- stackApply(r.feat.reg,c(1,1,1),fun=reg.id) #dif block dif ids (1 layer of id)

#Convert into data frame
coord <- round(coordinates(sample.unit.id),2)
sample.id <- getValues(sample.unit.id)
sample.area <- getValues(raster.ecoReg.area.l)
sample.unit.df <- as.data.frame(cbind(coord,sample.id,sample.area))
sample.unit.df <- na.omit(sample.unit.df)
names(sample.unit.df)
plot(y~x, data=sample.unit.df,cex=0.1)

##Create sample units of only a certain size
max.size <- 250000

sample.id <- function(x){if(x["sample.area"] < max.size){strsplit(as.character(x["sample.id"]),".",fixed=TRUE)[[1]][1]}else{x["sample.id"]}} #if ecoregion larger-split it up
blockingsampleunits <- apply(sample.unit.df,1,sample.id)
blockingsampleunits <- cbind(sample.unit.df,blockingsampleunits)
blockingsampleunits <- blockingsampleunits[,c("x","y","blockingsampleunits")]
colnames(blockingsampleunits) <- c("lon","lat","id.sample")

# Save the sample units for later use (Creating the units everytime takes too long)
write.csv(blockingsampleunits, "data/blockingsampleunits.csv", row.names=F)
library(ggplot2)
blockingsampleunits$id.sample <- as.numeric(blockingsampleunits$id.sample)
ggplot() + geom_raster(data=blockingsampleunits, aes(x=lon,y=lat, fill=id.sample))
