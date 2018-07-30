# Example for Jags-Yord-Xnom2grp-MrobustHet.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load The data file 
#....................................................................
myDataFrame = read.csv( file="OrdinalProbitData1.csv" )
yName="Y"
xName="X"
fileNameRoot="OrdinalProbitData1-"
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Yord-Xnom2grp-MnormalHet.R")
#------------------------------------------------------------------------------- 
# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=12000 , thinSteps=3 , 
#                    numSavedSteps=5000 , thinSteps=1 , 
                    saveName=fileNameRoot ) 
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
# Get all parameter names:
parameterNames = varnames(mcmcCoda) 
# for ( parName in parameterNames ) {
#   diagMCMC( codaObject=mcmcCoda , parName=parName , 
#                 saveName=fileNameRoot , saveType=graphFileType )
# }
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=c(-0.2,0.2) , 
                        RopeSdDiff=c(-0.2,0.2) , RopeEff=c(-0.1,0.1) , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
          RopeMuDiff=c(-0.2,0.2) , RopeSdDiff=c(-0.2,0.2), RopeEff=c(-0.1,0.1) , 
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

# NHST tests:
# t test difference of means:
tInfo = t.test(myDataFrame[,yName]~myDataFrame[,xName])
show(tInfo)
# effect size and confidence interval on effect size:
library(compute.es)
tesInfo = tes( t=abs(tInfo$stat) ,
               n.1=sum(as.numeric(myDataFrame[,xName])==1) ,
               n.2=sum(as.numeric(myDataFrame[,xName])==2) , dig=4 )
cat("\n\n\n")
show( tesInfo[c("pval.d","l.d","d","u.d")] )
cat("\n\n")
# variances and test of ratio of variances:
show(aggregate(myDataFrame[,yName]~myDataFrame[,xName],data=myDataFrame,sd))
show(var.test(myDataFrame[,yName]~myDataFrame[,xName],data=myDataFrame))
