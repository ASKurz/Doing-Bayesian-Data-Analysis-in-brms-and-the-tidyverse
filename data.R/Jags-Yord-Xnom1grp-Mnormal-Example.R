# Example for Jags-Yord-Xnom1grp-Mnormal.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
myDataFrame = read.csv( file="OrdinalProbitData-1grp-1.csv" )
yName="Y"
compVal=2.0
fileNameRoot="OrdinalProbitData-1grp-1-"
#.....................................................................
# myDataFrame = read.csv( file="OrdinalProbitData-1grp-2.csv" )
# yName="Y"
# compVal=4.0
# fileNameRoot="OrdinalProbitData-1grp-2-"
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Yord-Xnom1grp-Mnormal.R")
#------------------------------------------------------------------------------- 
# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , 
                    numSavedSteps=15000 , thinSteps=10 , 
#                    numSavedSteps=5000 , thinSteps=2 , 
                    saveName=fileNameRoot ) 
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
# Get all parameter names:
parameterNames = varnames(mcmcCoda) 
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
                saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , compVal=compVal , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , 
          compVal=compVal , 
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

# NHST tests:
# t test of mean:
tInfo = t.test(myDataFrame[,yName],mu=compVal)
show(tInfo)
tEffSz = (mean(myDataFrame[,yName])-compVal)/sd(myDataFrame[,yName])
show(tEffSz)
show(sd(myDataFrame[,yName]))

