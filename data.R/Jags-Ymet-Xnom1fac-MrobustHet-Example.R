# Example for Jags-Ymet-Xnom1fac-MrobustHet.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load The data file 

# Comment out or uncomment one of the data sections below. 
# In RStudio, select section then press Cntrl-Shift-C

# myDataFrame = read.csv( file="FruitflyDataReduced.csv" )
# # Specify the column names in the data file relevant to the analysis:
# yName="Longevity" 
# xName="CompanionNumber" 
# # Specify desired contrasts.
# # Each main-effect contrast is a list of 2 vectors of level names, 
# # a comparison value (typically 0.0), and a ROPE (which could be NULL):
# contrasts = list( 
#   list( c("Pregnant1","Pregnant8") , c("None0") , compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c("Pregnant1","Pregnant8","None0") , c("Virgin1","Virgin8") , 
#         compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c("Pregnant1","Pregnant8","None0") , c("Virgin1") , 
#         compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c("Virgin1") , c("Virgin8") , compVal=0.0 , ROPE=c(-1,1) ) 
# )
# # Specify filename root and graphical format for saving output.
# # Otherwise specify as NULL or leave saveName and saveType arguments 
# # out of function calls.
# fileNameRoot = "FruitflyData-RobustHet-" 
# graphFileType = "eps" 


myDataFrame = read.csv( file="NonhomogVarData.csv" )
yName="Y" 
xName="Group" 
contrasts = list( 
  list( c("D") , c("A") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("C") , c("B") , compVal=0.0 , ROPE=c(-1,1) ) 
)
fileNameRoot = "NonhomogVarData-RobustHet-" 
graphFileType = "eps" 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom1fac-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=11000 , thinSteps=20 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("aSigma","b0","b[1]","ySigma[1]", "nu","ySigmaMode","ySigmaSD") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , xName=xName ,
                        contrasts=contrasts , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , xName=xName ,
          contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
