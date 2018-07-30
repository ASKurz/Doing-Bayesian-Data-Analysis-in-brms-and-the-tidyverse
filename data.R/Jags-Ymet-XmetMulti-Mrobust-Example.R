# Example for Jags-Ymet-XmetMulti-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
# UNCOMMENT ONE OF THE FOLLOWING SECTIONS (In RStudio, select and ctrl-shift-C)
#.............................................................................
# # Two predictors:
myData = read.csv( file="Guber1999data.csv" )
yName = "SATT" ; xName = c("Spend","PrcntTake")
fileNameRoot = "Guber1999data-Jags-" 
numSavedSteps=15000 ; thinSteps=5
#.............................................................................
# # Single predictor:
# myData = read.csv( file="Guber1999data.csv" )
# yName = "SATT" ; xName = c("Spend")
# fileNameRoot = "Guber1999data-SinglePred-" 
# numSavedSteps=15000 ; thinSteps=5
#.............................................................................
# # Two predictors with redundant predictor:
# myData = read.csv( file="Guber1999data.csv" )
# PropNotTake = (100-myData[,"PrcntTake"])/100
# myData = cbind( myData , PropNotTake )
# yName = "SATT" ; xName = c("Spend","PrcntTake","PropNotTake")
# fileNameRoot = "Guber1999data-Jags-Redund-" 
# numSavedSteps=15000 ; thinSteps=15
#.............................................................................
# # Two predictors with two redundant predictors:
# myData = read.csv( file="Guber1999data.csv" )
# PropNotTake = (100-myData[,"PrcntTake"])/100
# Partic = myData[,"PrcntTake"]/10 
# myData = cbind( myData , PropNotTake , Partic )
# yName = "SATT" ; xName = c("Spend","PrcntTake","PropNotTake","Partic")
# fileNameRoot = "Guber1999data-Jags-Redund2-" 
# numSavedSteps=15000 ; thinSteps=15
#.............................................................................
# # Four predictors:
# myData = read.csv( file="Guber1999data.csv" )
# yName = "SATT" ; xName = c("Spend","PrcntTake","StuTeaRat","Salary")
# fileNameRoot = "Guber1999data-Jags-4X-" 
# numSavedSteps=15000 ; thinSteps=20 
#.............................................................................
# # Append columns of random predictors:
# myData = read.csv( file="Guber1999data.csv" )
# standardize = function(v){(v-mean(v))/sd(v)}
# Ny=nrow(myData)
# NxRand = 12
# set.seed(47405)
# for ( xIdx in 1:NxRand ) {
#   xRand = standardize(rnorm(Ny))
#   myData = cbind( myData , xRand )
#   colnames(myData)[ncol(myData)] = paste0("xRand",xIdx)
# }
# yName = "SATT" ; xName = c("Spend","PrcntTake", paste0("xRand",1:NxRand) )
# fileNameRoot = "Guber1999data-Jags-RandX-" 
# numSavedSteps=15000 ; thinSteps=5
#.............................................................................
# # Two predictors with interaction.
# # ** For two predictors with interaction, also uncomment graph command at end.
# # Append named column with interaction product:
# myData = cbind( myData , SpendXPrcnt=myData[,"Spend"]*myData[,"PrcntTake"] )
# yName = "SATT" ; xName = c("Spend","PrcntTake","SpendXPrcnt")
# fileNameRoot = "Guber1999data-Jags-Inter-" 
# numSavedSteps=15000 ; thinSteps=50
#.............................................................................
# # Introductory figures of chapter:
# myData = read.csv( file="MultLinRegrPlotUnif.csv" )
# myData = myData[101:150,]
# yName = "y" ; xName = c("x1","x2")
# fileNameRoot = "MultLinRegrPlotUnif-" 
# numSavedSteps=11000 ; thinSteps=2
#.............................................................................
# myData = read.csv( file="MultLinRegrPlotUnif.csv" )
# myData = myData[101:150,]
# yName = "y" ; xName = c("x1")
# fileNameRoot = "MultLinRegrPlotUnif-x1only-" 
# numSavedSteps=11000 ; thinSteps=2
#.............................................................................
# myData = read.csv( file="MultLinRegrPlotCorr.csv" )
# yName = "y" ; xName = c("x1","x2")
# fileNameRoot = "MultLinRegrPlotCorr-" 
# numSavedSteps=11000 ; thinSteps=2
#.............................................................................
# myData = read.csv( file="MultLinRegrPlotCorr.csv" )
# yName = "y" ; xName = c("x1")
# fileNameRoot = "MultLinRegrPlotCorr-x1only-" 
# numSavedSteps=11000 ; thinSteps=2
#.............................................................................

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-XmetMulti-Mrobust.R")
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
