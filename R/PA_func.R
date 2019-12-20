#############################################################################################
#             Select pseudo absences based on the distance to a species' range              #
#                 Using Naiara's distance weighting "One over distance 2"                   #
#                                      July 2014                                            #
#############################################################################################


#-#-# Function for absence selection #-#-#
PA.calc <- function (x){
  filename <- paste(strsplit(x,".",fixed=TRUE)[[1]][1])
  if(!file.exists(paste0(filetest,filename,"_PA1.Rdata"))){
    spData_alldistance <- get(load(paste0(spDistDir,x)))
    spDistr <- raster(paste0(spPresDir, filename,".tif"))
    coord <- round(coordinates(spDistr),4)
    presence <- getValues(spDistr)
    spDistr <- (as.data.frame(cbind(coord,presence)))
    spDistr[is.na(spDistr)] <- 0
    NP <- nrow(spDistr[which(spDistr$presence==1),])  # number of presences
    
    if((NP >= 10) == TRUE) { #Skip species with less than 10 presences
      
      if((NP > 35000) == FALSE){ #Return names for very widespread species where presence equals absence wont work
        
        speciesData <- merge(spData_alldistance,spDistr)
        allCoor <-  speciesData[,c("x","y")]
        distData <- speciesData[which(speciesData$presence == 0),] # Select the absence rows (contain distance info)
        
        if((NP < 1001) == TRUE) { #Depending which species need to be selected - use 1000 absences or NP!
          
          allData <- lapply(seq(1,10,1), function(i,NP,speciesData){
            smpl_distance <- distData[sample.int(nrow(distData),1000, prob=distData$OneOverDist2,  replace=F),c("x","y")] # randomly pick absence cell weighted by 1/D
            smpl_distance$presence <- 0
            smpl_distance <- rbind(smpl_distance,speciesData[which(speciesData$presence == 1),c("x","y","presence")])
            smpl_distance$sample_id <- i
            return(smpl_distance)
          },NP=NP,speciesData=speciesData)
          names(allData) <- c("PA1", "PA2", "PA3", "PA4", "PA5", 
                              "PA6", "PA7", "PA8", "PA9", "PA10")
          save(allData,file=paste0(filetest,filename,"_PA.Rdata"), compress="xz")   
          
        } else {
          allData <- lapply(seq(1,10,1), function(i,NP,speciesData){
            smpl_distance <- distData[sample.int(nrow(distData),NP, prob=distData$OneOverDist2,  replace=T),c("x","y")] # randomly pick absence cell weighted by 1/D
            smpl_distance$presence <- 0
            smpl_distance <- rbind(smpl_distance,speciesData[which(speciesData$presence == 1),c("x","y","presence")])
            smpl_distance$sample_id <- i
            return(smpl_distance)
          },NP=NP,speciesData=speciesData)
          names(allData) <- c("PA1", "PA2", "PA3", "PA4", "PA5", 
                              "PA6", "PA7", "PA8", "PA9", "PA10")
          save(allData,file=paste0(filetest, filename,"_PA.Rdata"), compress="xz")   
        }
      } else {
        wide.range <- c(wide.range, filename)
        print(filename) 
      }
    } else {
      print("no presence")
    }  
  } else {
    print("done")
  }
}