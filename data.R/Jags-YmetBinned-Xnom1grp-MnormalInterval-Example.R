# Example for Jags-YmetBinned-Xnom1grp-MnormalInterval.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

#.............................................................................
# Create random data for illustration:
set.seed(47401)
y = rnorm( 500 )
trueM = 100
trueSD = 15
y = (y-mean(y))/sd(y) * trueSD + trueM
# Corresponding bin values:
binThresholds = c( trueM-0.4*trueSD , trueM , trueM+0.4*trueSD )
ybin = y # placeholder
for ( i in 1:length(y) ) {
  ybin[i] = which( y[i] > c(-Inf,binThresholds) & y[i] < c(binThresholds,Inf)  , 
                   arr.ind=TRUE ) - 1
}
# Now change y values in some bins to NA:
yOrig = y # just for illustration
y[ ybin==1 ] = NA
y[ ybin==3 ] = NA
# Put y, ybin, and thresholds into a data frame:
myData = data.frame( yOrig=yOrig , ybin=ybin , y=y ,
                     matrix( rep( binThresholds , length(y) ) , 
                             byrow=TRUE , nrow=length(y) ,
                             dimnames = list( NULL ,
                                              paste0("thresh",1:length(binThresholds))
                                              ) ) )
write.csv( myData , file="CensoredNormalData.csv", row.names=FALSE )
rm(list=ls())  # Careful! This clears all of R's memory!
#.............................................................................

myData = read.csv("CensoredNormalData.csv")
yName="y" 
yBinName="ybin" 
threshVecName=c("thresh1","thresh2","thresh3")
fileNameRoot = "CensoredNormalData-ImputeData-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-YmetBinned-Xnom1grp-MnormalInterval.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( data=myData , yName=yName , yBinName=yBinName , 
                    threshVecName=threshVecName ,
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
summaryInfo = smryMCMC( mcmcCoda , 
                        compValMu=100.0 , #ropeMu=c(99.0,101.0) ,
                        compValSigma=15.0 , #ropeSigma=c(14,16) ,
                        compValEff=0.0 , #ropeEff=c(-0.1,0.1) ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , yName=yName , yBinName=yBinName , 
          threshVecName=threshVecName , 
          compValMu=100.0 , #ropeMu=c(99.0,101.0) ,
          compValSigma=15.0 , #ropeSigma=c(14,16) ,
          compValEff=0.0 , #ropeEff=c(-0.1,0.1) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
#------------------------------------------------------------------------------- 

# Now do again but omitting censored data:
rm(list=ls())  # Careful! This clears all of R's memory!
myData = read.csv("CensoredNormalData.csv")
myData=myData[,"y"] 
myData=myData[is.finite(myData)]
fileNameRoot = "CensoredNormalData-OmitCensor-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom1grp-Mnormal.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( data=myData , numSavedSteps=20000 , saveName=fileNameRoot )
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
                        compValMu=100.0 , #ropeMu=c(99.0,101.0) ,
                        compValSigma=15.0 , #ropeSigma=c(14,16) ,
                        compValEff=0.0 , #ropeEff=c(-0.1,0.1) ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , 
          compValMu=100.0 , #ropeMu=c(99.0,101.0) ,
          compValSigma=15.0 , #ropeSigma=c(14,16) ,
          compValEff=0.0 , #ropeEff=c(-0.1,0.1) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 


