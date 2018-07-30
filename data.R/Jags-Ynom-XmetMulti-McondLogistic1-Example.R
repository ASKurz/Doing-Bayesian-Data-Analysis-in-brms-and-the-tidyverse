#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#.............................................................................
myData = read.csv( file="CondLogistRegData1.csv" )
yName = "Y" ; xName = c("X1","X2")
fileNameRoot = "CondLogistRegData1-Model1-" 
numSavedSteps=5000 ; thinSteps=5
#.............................................................................
# myData = read.csv( file="CondLogistRegData2.csv" )
# yName = "Y" ; xName = c("X1","X2")
# fileNameRoot = "CondLogistRegData2-Model1-" 
# numSavedSteps=5000 ; thinSteps=5
#.............................................................................
# myData = read.csv( file="SoftmaxRegData1.csv" )
# yName = "Y" ; xName = c("X1","X2")
# fileNameRoot = "SoftmaxRegData1-Model1-" 
# numSavedSteps=5000 ; thinSteps=5
#.............................................................................
# myData = read.csv( file="SoftmaxRegData2.csv" )
# yName = "Y" ; xName = c("X1","X2")
# fileNameRoot = "SoftmaxRegData2-Model1-" 
# numSavedSteps=5000 ; thinSteps=5
#.............................................................................
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ynom-XmetMulti-McondLogistic1.R")
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
for ( parName in grep( "^beta" , parameterNames , value=TRUE ) ) {
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
