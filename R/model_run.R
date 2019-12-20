model_run <- function(sp){
  spname <- basename(sp)
  print(spname)
  clim.var <- as.character(unlist(climCombs))
  climVarName <- paste(climCombs,collapse="_")
  name <- unlist(lapply(sp, function(x) strsplit(basename(x),split=".",fixed=T)[[1]][1]))
  
  ## Get species data
  spdata <- get(load(paste0(sourceObs,"/", spname, ".Rdata")))
  
  # Select the presence cells again and count them 
  (ncells.pres <- nrow(spdata[[1]][spdata[[1]]$presence==1,]))
  
  if(ncells.pres >= 50){   ## Skip restricted range species
    # Run model for each PA set
    mod <- lapply(1:10, function(y){
      PA <- paste0(name, y) 
      species.data <- spdata[[y]][,c("x","y","presence")]
      species.data <- na.omit(species.data)
      spPseudoRep <- merge(species.data,climVar,by=c("x","y"),all.x=T)
      spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, 
                                   dplyr::matches("block"), dplyr::one_of(clim.var))
      spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
      if(model_type == "RF"){spPseudoRep <- na.omit(spPseudoRep)}
      
      ## Sum presences per block
      block.sum.all <- aggregate(presence~block,data=spPseudoRep,FUN=function(x)length(x)) # Overall length of each block (presences + absences)
      block.sum.pres <- aggregate(presence~block,data=spPseudoRep,FUN=sum) # Sum presences per block 
      block.sum <- merge(block.sum.all,block.sum.pres,by="block",all.x=T)
      
      block.not.zero <- block.sum[block.sum[,3] > 0 & block.sum[,2] >= 10,] # Define how many presences and absences each block must have to be considered
      num.block.not.zero <- nrow(block.not.zero) # Number of blocks that fulfill the presence absence criteria
      
      ## Remove blocks that have zero presences
      if(num.block.not.zero > 1){  
        # There must be at least one block with presences
        
        # Include only the blocks that have presences
        block.include <- block.not.zero[,1]
        
        ## Model function GAM
        if(model_type == "GAM"){
          GAM_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, PA=PA, 
                  clim.var=clim.var, fx=FALSE, k=-1, bs="tp",blocks=block.include)
        } else if(model_type == "GBM"){
          ## Model function GBM
          GBM_eco(data.model=spPseudoRep, outDir=resultsPath, plotPath=plotPath,
                  species=sp, PA=PA, clim.var=clim.var, blocks=block.include)
        } else if(model_type == "RF"){
          ## Model function RandomForest
          RF_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                 clim.var=clim.var, blocks=block.include, PA = PA)
        } else if(model_type == "MaxEnt"){
          ## Model function RandomForest
          MaxEnt_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                     clim.var=clim.var, blocks=block.include, PA = PA)
        }
      } else{
        mod <- NULL
      }
    })
    mod <- Filter(Negate(is.null), mod)
    save(mod, file=paste(resultsPath, "/", sp,"_", paste(clim.var,collapse="_"), "_model_output_", 
                           model_type, "_Eco_block.RData",sep=""), compress="xz")
  } else if(ncells.pres >= 10){
    # Run model for each PA set
    mod <- lapply(1:10, function(y){
      PA <- paste0(name, y) 
      species.data <- spdata[[y]][,c("x","y","presence")]
      species.data <- na.omit(species.data)
      spPseudoRep <- merge(species.data,climVar,by=c("x","y"),all.x=T)
      spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, dplyr::matches("block"), dplyr::one_of(clim.var))
      spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
      if(model_type == "RF"){spPseudoRep <- na.omit(spPseudoRep)}
      
      if(model_type == "GAM"){
        # Model function GAM
        GAM_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                  PA=PA, clim.var=clim.var, fx=FALSE, k=-1, bs="tp")
      } else if(model_type == "GBM"){
        GBM_split(data.model=spPseudoRep, outDir=resultsPath, plotPath=plotPath, 
                  species=sp, PA=PA, clim.var=clim.var)
      } else if(model_type == "RF"){
        RF_split(data.model=spPseudoRep, 
                 outDir=resultsPath, 
                 species=sp, PA=PA, clim.var=clim.var)
      } else if(model_type == "MaxEnt"){
        MaxEnt_split(data.model=spPseudoRep, outDir=resultsPath, 
                     species=sp, PA=PA, clim.var=clim.var)
      }
    })
    mod <- Filter(Negate(is.null), mod)
    if(length(mod) > 5){
    save(mod, file=paste(resultsPath, "/", sp,"_", paste0(clim.var, collapse="_"), 
                         "_model_output_", model_type, "_30_70.RData", sep=""), compress="xz")
    }
  }
  if(ncells.pres >= 10){
    if(length(mod) <= 5){
      #### Run 30-70 models for species where Ecoblocking did not create more than 5 models
      if(ncells.pres >= 50){
        mod <- lapply(1:10, function(y){
          PA <- paste0(name, y) 
          species.data <- spdata[[y]][,c("x","y","presence")]
          species.data <- na.omit(species.data)
          spPseudoRep <- merge(species.data,climVar,by=c("x","y"),all.x=T)
          spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, dplyr::matches("block"), dplyr::one_of(clim.var))
          spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
          if(model_type == "RF"){spPseudoRep <- na.omit(spPseudoRep)}
          
          if(model_type == "GAM"){
            # Model function GAM
            GAM_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, PA=PA, 
                      clim.var=clim.var, fx=FALSE, k=-1, bs="tp")
          } else if(model_type == "GBM"){
            GBM_split(data.model=spPseudoRep, outDir=resultsPath, plotPath=plotPath,
                      species=sp, PA=PA, clim.var=clim.var)
          } else if(model_type== "RF"){
            RF_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                     PA=PA, clim.var=clim.var)
          } else if(model_type == "MaxEnt"){
            MaxEnt_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                         PA=PA, clim.var=clim.var)
          }
        })
        mod <- Filter(Negate(is.null), mod)
        if(length(mod) > 5){
        save(mod, file=paste(resultsPath, "/", sp,"_", paste0(clim.var, collapse="_"), 
                             "_model_output_", model_type, "_30_70_MissingEco.RData",
                             sep=""), compress="xz")
        }
      }
    }
  }
  gc()
  #raster::removeTmpFiles(h=0.1)
}
