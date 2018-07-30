
source("minNforHDIpower.R")

# Specify desired powers:
desPow = seq(.7,.9,by=0.1)
# Specify generating modes:
genMod = seq(0.60,0.85,by=0.05)

# HDI excludes ROPE around nullVal:
# Declare matrix for storing results:
sampSizeMatrix = matrix( NA , nrow=length(desPow) , ncol=length(genMod) ,
                         dimnames=list("Power"=desPow,"Gen.Mode"=genMod) )
# Compute sample size for all combinations of desired powers and modes:
for ( desPowIdx in 1:length(desPow) ) {
  for ( genModIdx in 1:length(genMod)  ) {
    sampSize = minNforHDIpower( genPriorMode=genMod[genModIdx], genPriorN=2000,
                                HDImaxwid=NULL, nullVal=0.5, ROPE=c(0.48,0.52),
                                desiredPower=desPow[desPowIdx] , 
                                audPriorMode=0.5 , audPriorN=2 ,
                                HDImass=0.95 , initSampSize=5 , verbose=FALSE )
    sampSizeMatrix[desPowIdx,genModIdx ] = sampSize
    show(sampSizeMatrix)
  }
}
sampSizeExclNull = sampSizeMatrix
show(sampSizeExclNull)

# HDI has maximum width:
# Declare matrix for storing results:
sampSizeMatrix = matrix( NA , nrow=length(desPow) , ncol=length(genMod) ,
                         dimnames=list("Power"=desPow,"Gen.Mode"=genMod) )
# Compute sample size for all combinations of desired powers and modes:
for ( desPowIdx in 1:length(desPow) ) {
  for ( genModIdx in 1:length(genMod)  ) {
    sampSize = minNforHDIpower( genPriorMode=genMod[genModIdx], genPriorN=10,
                                HDImaxwid=0.20, nullVal=NULL, ROPE=NULL,
                                desiredPower=desPow[desPowIdx] , 
                                audPriorMode=0.5 , audPriorN=2 ,
                                HDImass=0.95 , initSampSize=50 , verbose=FALSE )
    sampSizeMatrix[desPowIdx,genModIdx ] = sampSize
    show(sampSizeMatrix)
  }
}
sampSizeMaxHdiWid = sampSizeMatrix
show(sampSizeExclNull)
show(sampSizeMaxHdiWid)



