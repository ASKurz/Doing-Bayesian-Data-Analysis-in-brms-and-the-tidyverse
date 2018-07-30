# Example for Stan-Ymet-Xnom2grp-MrobustHet.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load The data file 
myDataFrame = read.csv( file="TwoGroupIQ.csv" )
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Stan-Ymet-Xnom2grp-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "TwoGroupIQrobustHet-Stan-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName="Score" , xName="Group" ,
                    numSavedSteps=20000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=c(-0.5,0.5) , 
                        RopeSdDiff=c(-0.5,0.5) , RopeEff=c(-0.1,0.1) , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName="Score" , xName="Group" , 
          RopeMuDiff=c(-0.5,0.5) , RopeSdDiff=c(-0.5,0.5), RopeEff=c(-0.1,0.1) , 
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
