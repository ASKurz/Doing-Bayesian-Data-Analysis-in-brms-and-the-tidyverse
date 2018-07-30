# Example for Jags-Ymet-XnomSplitPlot-MnormalHom.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#Load The data file 

fileNameRoot = "SplitPlotAgriData-" 
graphFileType = "eps" 
myDataFrame = read.csv( "SplitPlotAgriData.csv" )
# Specify the column names in the data file relevant to the analysis:
yName="Yield" 
xBetweenName="Till" 
xWithinName="Fert" 
xSubjectName="Field"
xBetweenContrasts = list( 
  list( c("Moldbrd") , c("Ridge") , compVal=0.0 , ROPE=c(-5,5) ) ,
  list( c("Moldbrd","Ridge") , c("Chisel") , compVal=0.0 , ROPE=c(-5,5) ) 
)
xWithinContrasts = list( 
  list( c("Deep","Surface") , c("Broad") , compVal=0.0 , ROPE=c(-5,5) ) 
)
xBetweenWithinContrasts = list(
  list( list( c("Chisel","Moldbrd") , c("Ridge") ) ,
        list(  c("Broad") , c("Deep","Surface") ) ,
        compVal=0.0 , ROPE=c(-5,5) ) 
)

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-XnomSplitPlot-MnormalHom.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , xBetweenName=xBetweenName , 
                    xWithinName=xWithinName , xSubjectName=xSubjectName ,
                    numSavedSteps=12000 , thinSteps=20 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("mSxW[1,1]","b0","bB[1]","bW[1]","bS[1]","bBxW[1,1]",
                   "sigma","sigmaB","sigmaW","sigmaS","sigmaBxW") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , xBetweenName=xBetweenName , 
                        xWithinName=xWithinName , xSubjectName=xSubjectName ,
                        xBetweenContrasts=xBetweenContrasts , 
                        xWithinContrasts=xWithinContrasts , 
                        xBetweenWithinContrasts=xBetweenWithinContrasts ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , xBetweenName=xBetweenName , 
          xWithinName=xWithinName , xSubjectName=xSubjectName ,
          xBetweenContrasts=xBetweenContrasts , 
          xWithinContrasts=xWithinContrasts , 
          xBetweenWithinContrasts=xBetweenWithinContrasts ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
