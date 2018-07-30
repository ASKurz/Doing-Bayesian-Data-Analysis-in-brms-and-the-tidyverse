# Example for Stan-Ymet-XmetSsubj-MrobustHierQuadWt.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):

# myData = read.csv( file="IncomeFamszState3yr.csv" , comment.char="#")
# xName = "FamilySize" ; yName = "MedianIncome" ; sName="State" ; wName="SampErr"
# fileNameRoot = "IncomeFamszState3yr-Stan-" 

 myData = read.csv( file="HierLinRegressData.csv" )
 xName = "X" ; yName = "Y" ; sName="Subj" ; wName=NULL
 fileNameRoot = "HierLinRegressData-Quad-Stan-" 

# myData = read.csv( file="BugsRatsData.csv" )
# xName = "Day" ; yName = "Weight" ; sName="Subj" ; wName=NULL
# fileNameRoot = "BugsRatsData-Quad-Stan-" 

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Stan-Ymet-XmetSsubj-MrobustHierQuadWt.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , 
                    xName=xName , yName=yName , sName=sName , wName=wName ,
                    numSavedSteps=12000 , thinSteps=5 , saveName=fileNameRoot )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
# #------------------------------------------------------------------------------- 
# # Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in c("beta0mu","beta1mu","beta2mu","nu","sigma",
                   "beta0[1]","beta1[1]","beta2[1]") ) {
 diagMCMC( codaObject=mcmcCoda , parName=parName , 
           saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , 
          xName=xName , yName=yName , sName=sName , wName=wName ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
