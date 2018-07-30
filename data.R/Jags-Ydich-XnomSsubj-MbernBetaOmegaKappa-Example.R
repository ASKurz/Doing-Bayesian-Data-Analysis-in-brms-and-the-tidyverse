# Example for Jags-Ydich-XnomSsubj-MbernBetaOmegaKappa.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Read The data file:
myData = read.csv("TherapeuticTouchData.csv")
yName = "y" # column name for 0,1 values
sName = "s" # column name for subject ID
# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "Jags-Ydich-XnomSsubj-MbernBetaOmegaKappa-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# # Read The data file:
# myData = read.csv("StormTressoldiDiRisio2010data.csv")
# yName = "Correct" # column name for 0,1 values
# sName = "Study" # column name for "subject" ID
# # Optional: Specify filename root and graphical format for saving output.
# # Otherwise specify as NULL or leave saveName and saveType arguments 
# # out of function calls.
# fileNameRoot = "StormTressoldiDiRisio2010-" 
# graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ydich-XnomSsubj-MbernBetaOmegaKappa.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( data=myData , sName=sName , yName=yName , 
                    numSavedSteps=20000 , saveName=fileNameRoot , thinSteps=10 )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names for reference
for ( parName in parameterNames[c(1:3,length(parameterNames))] ) { 
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
                saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , compVal=0.5 , 
                        diffIdVec=c(1,14,28),  # Therapeutic touch
                        # diffIdVec=c(38,60,2),  # ESP Tressoldi et al.
                        compValDiff=0.0 ,
                        saveName=fileNameRoot )
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , sName=sName , yName=yName , 
          compVal=0.5 , #rope=c(0.45,0.55) , # Therapeutic touch
          diffIdVec=c(1,14,28),              # Therapeutic touch
          # compVal=0.25 , #rope=c(0.22,0.28) , # ESP Tressoldi et al.
          # diffIdVec=c(38,60,2),               # ESP Tressoldi et al.
          compValDiff=0.0, #ropeDiff = c(-0.05,0.05) ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
