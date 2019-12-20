#' ## Blocking the global climate data by ecoregion 
#' Based on Robies' blocking method GCB paper 2013 
#' July 2014

#' See create_blockingsampleunits.R for creating the BlockingSampleUnits.csv file.

#' ## Create the blocks using the baseline climate data

#' Read in the data

# Read in the csv with the Blocking Samples
sample.units.id<-read.csv("extdata/BlockingSampleUnits.csv")[,2:4]
colnames(sample.units.id) <- c("x","y","id.sample")

# Read in the climate data
climData <- read.csv("extdata/bioclim_1995.csv")[,-1]
names(climData)
out<-climData
sample.units.id <- merge(sample.units.id,out,by=c("x","y"),all.x=FALSE)
names(sample.units.id)
head(sample.units.id)

#' Do the blocking

# Aggregate to sample regions 
# !!Chose your climate variables you want to model with and block with those!!
sample.unit.mean <- aggregate(cbind(bio1,bio4, bio12, bio15)~id.sample,data=sample.units.id,mean)
colnames(sample.unit.mean) <- c("id.sample","BIO1m","BIO4m","BIO12m","BIO15m")
head(sample.unit.mean)
sample.unit.var <- aggregate(cbind(bio1,bio4,bio12,bio15)~id.sample,data=sample.units.id,var)
colnames(sample.unit.var) <- c("id.sample","BIO1v","BIO4v","BIO12v","BIO15v")
head(sample.unit.var)
sample.unit.var[is.na(sample.unit.var)] <- 0
sample.unit.all <- merge(sample.unit.mean,sample.unit.var,by="id.sample",all.x=TRUE)
length(sample.unit.all[,1])
# names(sample.unit.all)
head(sample.unit.all)

#Create blocks, dividing polygons into orthogonal blocks on the basis of the climate data.
blocks <- block(sample.unit.all, n.tr=10, id.vars='id.sample') ## is this where you decide the number of blockds to use i.e. n.tr=5 or 10
blocks <- assignment(blocks, namesCol=as.character(1:10))$assg[[1]][1:10] ## assign subpols to one of 10 blocks 
blocks <- Reduce(rbind, mapply(function(id.sample, block) data.frame(id.sample, block), id.sample=blocks, block=as.list(1:10),SIMPLIFY=F) )## turn this into a dataframe 
blocks <- blocks[!is.na(blocks$id.sample),] ## remove any polygons that haven't been assigned to a block (we deal with this later).
head(blocks)
clidat <- merge(sample.unit.all, blocks, by=c('id.sample'), all=T) 
head(clidat)
clidat$block[is.na(clidat$block)] <- sample(1:10,1)
clidat <- clidat[,c(1,10)] # id sample and block ## this value will need to change depending on how many variables I have
names(clidat)

t <- merge(sample.units.id, clidat,by="id.sample",all.x=TRUE)
head(t)
levelplot(block~x+y,data=t)
names(sample.units.id)

write.csv(t, "extdata/Blocking_SU_1995.csv")