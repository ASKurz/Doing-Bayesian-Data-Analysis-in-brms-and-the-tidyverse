#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# #.............................................................................
# # Two predictors:
# myData = read.csv( file="HtWtData110.csv" )
# yName = "male" ; xName = c("weight","height")
# fileNameRoot = "HtWtData110-" 
# numSavedSteps=15000 ; thinSteps=2
# #.............................................................................
# Only one predictor:
myData = read.csv( file="HtWtData110.csv" )
yName = "male" ; xName = c("weight")
fileNameRoot = "HtWtData110-weightOnly-" 
numSavedSteps=15000 ; thinSteps=2
# #.............................................................................
# # Add some outliers:
# outlierMat = matrix( c(
#   190,74,0,
#   230,73,0,
#   120,59,1,
#   150,58,1 ) , ncol=3 , byrow=TRUE , 
#   dimnames= list( NULL , c("weight","height","male") ) )
# myData = rbind( myData , outlierMat )
#.............................................................................
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ydich-XmetMulti-Mlogistic.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                    saveName=fileNameRoot )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
