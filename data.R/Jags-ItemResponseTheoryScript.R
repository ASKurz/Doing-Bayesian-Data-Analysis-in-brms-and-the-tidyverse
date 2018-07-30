# Jags-ItemResponseTheoryScript.R 
# John Kruschke, November 2015.
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

# Load the functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.

fileNameRoot = "Jags-ItemResponseTheory-"

#want = c("") # desired packages
#have = want %in% rownames(installed.packages())
#if ( any(!have) ) { install.packages( want[!have] ) }

# Define logistic function for creating data and for later plotting
logistic = function( x , g=1 , t=0 ) {
  return( 1/(1+exp(-(g*(x-t)))) )
}

#======================================================================
# Create some data for IRT:
createData=TRUE
if ( createData ) {
  # Specify item and subject properties:
  Nitem = 12
  Nsubj = 120
  set.seed(47405)
  itemDifficulty = runif(Nitem,0,100)
  itemDiscrimination = runif(Nitem,0.02,0.2)
  subjAbility = rnorm(Nsubj,50,25)
  Nitem = length(itemDifficulty)
  Nsubj = length(subjAbility)
  # Frame for holding data:
  theData = data.frame( Subj=rep(0,Nsubj*Nitem) , 
                        Item=rep(0,Nsubj*Nitem) , 
                        Correct=rep(0,Nsubj*Nitem) )
  # Generate random 0,1 responses:
  rowIdx = 0
  for ( subjIdx in 1:Nsubj ) {
    for ( itemIdx in 1:Nitem ) {
      rowIdx = rowIdx+1
      pCorrect = logistic( x=subjAbility[subjIdx] , 
                           g=itemDiscrimination[itemIdx] ,
                           t=itemDifficulty[itemIdx] )
      theData[rowIdx,"Subj"] = subjIdx
      theData[rowIdx,"Item"] = itemIdx
      theData[rowIdx,"Correct"] = sample( c(0,1) , size=1 , 
                                          prob=c(1-pCorrect,pCorrect) )
    }
  }
  # Append some incompete subject data:
  subjIdx = subjIdx + 1
  rowIdx=rowIdx+1
  theData[rowIdx,c("Subj","Item","Correct")]=c(subjIdx,7,1)
  subjIdx = subjIdx + 1
  rowIdx=rowIdx+1
  theData[rowIdx,c("Subj","Item","Correct")]=c(subjIdx,7,0)
  # Make sure the subj and item ID are consecutive integers:
  theData$Subj = factor( theData$Subj )
  theData$Item = factor( theData$Item )
  theData$Correct = as.integer( theData$Correct )
  # Save the data so it can be read like a real data file:
  write.csv( theData , file="ItemResponseTheoryData.csv" , row.names=FALSE )
}
# End of creating data
#======================================================================

theData = read.csv( file="ItemResponseTheoryData.csv" )
# Make sure the subj and item ID are consecutive integers:
theData$Subj = factor( theData$Subj )
theData$Item = factor( theData$Item )
theData$Correct = as.integer( theData$Correct )
#str(theData)

#######################################################################
# THE CODE BELOW ASSUMES THAT 
# THE DATA COLUMN NAMES ARE "Subj" "Item" and "Correct"
#######################################################################

# Identify difficult and easy items for calibration:
subjPropCorr = aggregate( Correct ~ Subj , data=theData , FUN="mean" ) 
itemPropCorr = aggregate( Correct ~ Item , data=theData , FUN="mean" ) 
itemDiffMaxIdx = which.min( # index of hardest not-impossible item
  itemPropCorr$Correct[ itemPropCorr$Correct > 0 ] )
itemDiffMinIdx = which.max( # index of easiest not-trivial item
  itemPropCorr$Correct[ itemPropCorr$Correct < 1 ] )

# Assemble data into list for JAGS:
y = theData[,"Correct"] 
itemID = as.numeric(factor(theData[,"Item"])) 
subjID = as.numeric(factor(theData[,"Subj"])) 
Nitem = length(unique(itemID)) 
Nsubj = length(unique(subjID)) 
Ntotal = nrow(theData) 
itemDiff = rep(NA,Nitem)   # itemDiff is estimated by JAGS
itemDiff[itemDiffMinIdx]=0   # itemDiff of easy item is set to 0
itemDiff[itemDiffMaxIdx]=100   # itemDiff of hard item is set to 100
itemDiffEstIdx = which(is.na(itemDiff)) # indices of remaining items
dataList = list(
  y=y , itemID=itemID , subjID=subjID , Nitem=Nitem , Nsubj=Nsubj , 
  Ntotal=Ntotal , itemDiff=itemDiff , itemDiffEstIdx=itemDiffEstIdx 
)

# Define the model:
modelString = "
model {
  for ( rowIdx in 1:Ntotal ) {
    y[rowIdx] ~ dbern( pCorr[rowIdx] )
    pCorr[rowIdx] <- ilogit(  # ilogit is logistic in JAGS
                             itemDisc[itemID[rowIdx]] 
                              * ( subjAbil[subjID[rowIdx]] 
                                  - itemDiff[itemID[rowIdx]] ) )
  }
  for ( subjIdx in 1:Nsubj ) {
    subjAbil[subjIdx] ~ dnorm( muAbil , 1/sigmaAbil^2 )
  }
  for ( itemIdx in itemDiffEstIdx ) { # two items are anchored
    itemDiff[itemIdx] ~ dnorm( muDiff , 1/sigmaDiff^2 )
  }
  for ( itemIdx in 1:Nitem ) {  
    itemDisc[itemIdx] ~ dnorm( muDisc , 1/sigmaDisc^2 )
  }
  muDisc ~ dnorm(0.1,0.1^2)    # <- 0.1
  sigmaDisc ~ dunif(0.001,0.1) # <- 0.05
  muAbil     ~ dnorm(50,1/25^2)  # <- 50
  sigmaAbil  ~ dunif(1,50)       # <- 50
  muDiff <- 50    # ~ dnorm(50,1/50^2)
  sigmaDiff <- 50 # ~ dunif(0.01,100)
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )


# Run the chains:
require(runjags)
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" , 
                        monitor=c("subjAbil","muAbil","sigmaAbil",
                                  "itemDiff","muDiff","sigmaDiff",
                                  "itemDisc","muDisc","sigmaDisc" ) , 
                        data=dataList ,  
                        #inits=initsList , 
                        n.chains=3 ,
                        adapt=1000 ,
                        burnin=2000 , 
                        sample=ceiling(30000/3) , # less for demo only
                        thin=10 ,                 # less for demo only
                        summarise=FALSE ,
                        plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

save( codaSamples , file=paste0(fileNameRoot,"codaSamples.Rdata") )

# Convergence diagnostics:
parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in c(
  "subjAbil[1]","subjAbil[2]","muAbil","sigmaAbil",
  "itemDisc[1]","itemDisc[2]","muDisc","sigmaDisc",
  "itemDiff[1]","itemDiff[2]"
) ) {
  diagMCMC( codaObject=codaSamples , parName=parName ) 
}


# Examine the posterior distribution:
mcmcMat = as.matrix(codaSamples)
chainLength = nrow(mcmcMat)
chainSubSteps = ceiling(seq(1,chainLength,length=20))

mcmcMean = apply( mcmcMat , 2 , mean )
itemInfo = cbind( 
  round( mcmcMean[ grep( "itemDiff" , names(mcmcMean) ) ] , 1 ) ,
  round( mcmcMean[ grep( "itemDisc" , names(mcmcMean) ) ] , 3 ) )
subjInfo = cbind( 
  round( mcmcMean[ grep( "subjAbil" , names(mcmcMean) ) ] , 1 ) )

# Plot item results:
abilRange = range( subjInfo ) # range of mean abilities
discRange = range( mcmcMat[,grep("itemDisc",colnames(mcmcMat))] )
xComb = seq(abilRange[1],abilRange[2],length=301)
for ( itemIdx in 1:Nitem ) {
  openGraph()
  layout( matrix(c(1,1,2,3),nrow=2,byrow=TRUE) , heights=c(3,2) )
  par( mar=c(4,4,4,1) , mgp=c(2.1,0.7,0) )
  par( mar=c(4,4,4,1) , mgp=c(2.1,0.7,0) )
  plot( -100,-1, xlim=abilRange , ylim=c(0,1) , xlab="Difficulty (or Ability)" ,
        ylab="Probability Correct" , main=bquote("Item "*.(itemIdx)) ,
        cex.main=2 , cex.lab=1.5 )
  for ( stepIdx in chainSubSteps ) {
    lines( xComb , logistic( xComb , 
                             g=mcmcMat[stepIdx,paste0("itemDisc[",itemIdx,"]")] ,
                             t=mcmcMat[stepIdx,paste0("itemDiff[",itemIdx,"]")] ) ,
           col="skyblue" )
  }
  par( mar=c(4,1,1,1) , mgp=c(2.1,0.7,0) )
  plotPost( mcmcMat[,paste0("itemDiff[",itemIdx,"]")] , 
            xlab="Item Difficulty" , xlim=abilRange )
  plotPost( mcmcMat[,paste0("itemDisc[",itemIdx,"]")] , 
            xlab="Item Discrimination" , xlim=discRange )
}

# Plot subject info:

# Plot first 25 subject abilities:
openGraph()
layout( matrix(1:25,nrow=5,byrow=TRUE) )
par( mar=c(3.5,2.5,2.5,1.5) , mgp=c(2.1,0.7,0) )
for ( subjIdx in 1:min(c(25,Nsubj)) ) { 
  abilRange = range( subjInfo ) # range of mean abilities
  plotPost( mcmcMat[,paste0("subjAbil[",subjIdx,"]")] , 
            xlab="Ability" , xlim=abilRange , main=paste("Subj",subjIdx) )
}

# Plot selected individual subject abilities:
openGraph(height=7,width=7)
layout( matrix(1:4,nrow=2) )
par( mar=c(3.5,2.5,2.5,1.5) , mgp=c(2.1,0.7,0) )
for ( subjIdx in c(10,122,17,121) ) { 
  abilRange = range( subjInfo ) # range of mean abilities
  plotPost( mcmcMat[,paste0("subjAbil[",subjIdx,"]")] , 
            xlab="Ability" , xlim=abilRange , main=paste("Subj",subjIdx) )
}

# explore differences, correlations of parameters...

plotCompare = function( parName1 , parName2 , 
                        mcmcM=mcmcMat , NstepsToPlot=900 ) {
  stepsToPlot = ceiling(seq(1,nrow(mcmcM),length=NstepsToPlot))
  openGraph(width=4,height=7)
  layout( matrix(1:2,nrow=2) )
  par( mar=c(3.5,3.5,1.5,1.5) , mgp=c(2.1,0.7,0) )
  plot( mcmcM[stepsToPlot,parName2] , mcmcM[stepsToPlot,parName1] , 
        xlab=parName2 , ylab=parName1 , col="skyblue" )
  abline(0,1,lty="dashed")
  par( mar=c(3.5,1.5,2.5,1.5) , mgp=c(2.1,0.7,0) )
  plotPost( mcmcM[,parName1] - mcmcM[,parName2] , compVal=0.0 ,
            xlab="Difference" , main=paste(parName1,"-",parName2) )
}

#plotCompare( "itemDiff[12]" , "itemDiff[4]" )  
#plotCompare( "itemDiff[5]" , "itemDiff[4]" )  
#plotCompare( "subjAbil[9]" , "subjAbil[20]" )  
#plotCompare( "itemDisc[10]" , "itemDisc[3]" )  
plotCompare( "itemDiff[10]" , "itemDiff[1]" )  
plotCompare( "itemDisc[10]" , "itemDisc[12]" )  
plotCompare( "subjAbil[14]" , "subjAbil[11]" )  



