#############################################################################################
#      Calculate Vincenty distance for each absence cell to the nearest presence cell       #
#              Vincenty takes into account that the earth is a spheroid oblate              #
#                 Using Naiara's distance weighting "One over distance 2"                   #
#                                      July 2014                                            #
#############################################################################################

#' Function to calculate one over distance squared for all absence land cells for each species
distance.calc <- function (sp){
  # Extract the species name to use later when saving the distance data
  sp.name <- paste(strsplit(basename(sp),"[_.]",fixed=FALSE)[[1]][1:2],collapse="_")
  
  if(!file.exists(paste0(filetest,sp.name,".Rdata"))){
    
    #ObDist <- get(load(sp))[,c("x","y","presence")] #to process r data files
    ObDist <- raster(sp)
    coord <- round(coordinates(ObDist),4)
    presence <- getValues(ObDist)
    ObDist <- (as.data.frame(cbind(coord,presence)))
    ObDist[is.na(ObDist)] <- 0
    #ggplot()+geom_raster(data=ObDist,aes(x=x,y=y,fill=presence)) #plot to double check if distribution looks right
    
    if((sum(ObDist$presence) > 0) == TRUE) { # Select species with more than 30 cells (or other threshold)
      
      # Create a window with the raster extent
      w <- owin(c(-180,180), c(-90,90))
      
      # Select the absence cells in one data frame 
      abs.sub <- subset(ObDist, presence == 0)
      # Change the absence data into a point file
      abs.sub.final<-as.ppp(abs.sub,w)
      # Select the presence cells in another data frame
      pres.sub <- subset(ObDist, presence == 1)
      # Change the presence data into a point file too
      pres.sub.final<-as.ppp(pres.sub,w)
      # Identify the nearest present cell for each absence cell
      dist.cell.coord <- nncross(abs.sub.final, pres.sub.final,what = c("dist", "which"))   
      
      # Use the location of the nearest presence cell for each absence cell in the presence data to extract the coordinate 
      # Then merge the absence cell coordinate with the nearest presence cell coordinate
      abs.df<-list(abs.sub.final$x, abs.sub.final$y, dist.cell.coord)
      abs.df<-do.call(cbind,abs.df)
      colnames(abs.df)[1:2] <- c("x","y")
      
      pres.df<-list(pres.sub.final$x, pres.sub.final$y)
      pres.df<-as.data.frame(do.call(cbind,pres.df))
      colnames(pres.df)[1:2] <- c("x.p","y.p")
      pres.df$which <- rownames(pres.df)
      
      absNearPres <- merge(abs.df,pres.df,by="which")
      
      # Update the column names of the data frame that contains absence coordinates and nearest presence coordinates
      absNearPres <- absNearPres[,c("x","y","dist","which","x.p","y.p")]
      colnames(absNearPres) <- c("lon1","lat1","dist","which","lon2","lat2")
      absNearPres1<-as.data.frame(absNearPres)
      absNearPres<- absNearPres[,c("lat1","lon1","lat2","lon2")]
      absNearPres<-as.matrix(absNearPres)
      
      # Calculate the Vincenty distance for between each absence cell and the nearest presence cell 
      VincentyDist<- SDMTools::distance(absNearPres, bearing = FALSE)
      # Divide the distance by 1000 since it is in meters
      VincentyDist$distance<-VincentyDist$distance/1000
      # Update the columnnames of the Vincenty distance file
      colnames(VincentyDist)<-c("y","x","lat2","lon2","distance")
      # Select only the needed absence coordinates and their distances
      VincentyDist<-VincentyDist[,c("x","y","distance")]
      
      # The Vincenty distance file contains only the absence cells, so I need to add the presence coordinates and set the distance for those to 0
      pres.df["presence"] <- 0 
      pres.df <- pres.df[,c(1,2,4)]
      
      # Make the columnnames of the presence cells the same as from the absence cells
      colnames(pres.df)<-c("x","y","distance")
      
      # Put presence and absence cells into one dataframe
      FinalDist<-rbind(VincentyDist,pres.df)
      
      # Add another column in which Naiara's distance measure ("One over distance 2") is calculated 
      FinalDist["OneOverDist2"] <- 1/FinalDist$distance^2
      FinalDist <- merge(land,FinalDist,by=c("x","y"),all.x=T)
      FinalDist <- FinalDist[,c("x","y","OneOverDist2")]
      
      # Save the resulting distance dataframe as Rdata file for each species
      save(FinalDist, file=paste0(resultspath, sp.name,".Rdata"),compress="xz")
      
      removeTmpFiles(h=6)
      
    } else {
      print(sp.name)
    }
  } else {
    print(sp.name)
  } 
}