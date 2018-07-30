# In the style of Jags-ExampleScript.R, which accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

# Load the functions used below. Retrieve files from
# https://sites.google.com/site/doingbayesiandataanalysis/software-installation
source("DBDA2E-utilities.R") # Must be in R's current working directory.

# Data from Yusuf et al. (1985). In Table 5.4 of Bayesian Data Analysis,
# http://www.stat.columbia.edu/~gelman/book/data/meta.asc
# study control.deaths control.total treated.deaths treated.total
# 1     3 39          3 38 
# 2     14 116       7 114 
# 3     11 93         5 69 
# 4     127 1520  102 1533 
# 5     27 365      28 355 
# 6     6 52          4 59 
# 7     152 939     98 945 
# 8     48 471      60 632 
# 9     37 282      25 278 
# 10    188 1921  138 1916 
# 11    52 583      64 873 
# 12    47 266      45 263 
# 13    16 293       9 291 
# 14    45 883      57 858 
# 15    31 147      25 154 
# 16    38 213      33 207 
# 17    12 122      28 251 
# 18    6 154        8 151 
# 19    3 134        6 174 
# 20    40 218      32 209 
# 21    43 364      27 391 
# 22    39 674      22 680 
# # Also available as follows:
# install.packages("BayesDA")
# library("BayesDA")
# data(meta)
# colnames(meta) = c( "StudyID" , "ControlDeaths" , "ControlTotal" , 
#                     "TreatedDeaths" , "TreatedTotal" )
# write.csv( meta , row.names=FALSE , file="MetaAnalysisBetaBlocker.csv" )

theData = read.csv( "MetaAnalysisBetaBlocker.csv" )
fileNameRoot="Jags-MetaAnalysisBetaBlockers" # For output file names.

# Display a table of data:
ControlProp = theData$ControlDeaths / theData$ControlTotal 
TreatedProp = theData$TreatedDeaths / theData$TreatedTotal 
MeanProp = ( ControlProp + TreatedProp )/2 
DiffProp = TreatedProp - ControlProp 
multT = TreatedProp / ControlProp
show(
  cbind( theData , 
         ControlProp=round(ControlProp,3) , 
         multT=round(multT,3)
         # TreatedProp , 
         # MeanProp = round(MeanProp,3) , 
         # DiffProp = round(DiffProp,3) 
  )
)

# Package the data for JAGS:
dataList = list(
  nS = length(theData$StudyID) ,
  zC = theData$ControlDeaths ,
  nC = theData$ControlTotal , 
  zT = theData$TreatedDeaths ,
  nT = theData$TreatedTotal 
)

# Define the JAGS model:
modelString = "
model {
  for ( s in 1:nS ) {
    zC[s] ~ dbin( thetaC[s] , nC[s] )
    zT[s] ~ dbin( thetaT[s] , nT[s] )
    thetaC[s] ~ dbeta( moC*(kappaC-2)+1 , (1-moC)*(kappaC-2)+1 )
    thetaT[s] <- thetaC[s] * multT[s]
    multT[s] ~ dgamma( multTSh , multTRa )
  }
  multTSh <- 1 + multTMo * multTRa
  multTRa <- ( ( multTMo + sqrt( multTMo^2 + 4*multTSD^2 ) ) 
                / ( 2*multTSD^2 ) )
  multTMo ~ dgamma( 2.618 , 1.618 )  # mode=1 , sd=1
  multTSD ~ dgamma( 1.925 , 4.625 )  # mode=0.2 , sd=0.3
  # multTMo ~ dgamma( 1.925 ,  0.925 )  # mode=1 , sd=1.5
  # multTMo ~ dgamma( 1.640 ,  0.640 )  # mode=1 , sd=2
  # multTSD ~ dgamma( 2.618 , 16.180 )  # mode=0.1 , sd=0.1
  # multTSD ~ dgamma( 2.618 ,  8.090 )  # mode=0.2 , sd=0.2
  moC ~ dbeta(1,1)
  kappaC <- kappaCMinusTwo + 2
  kappaCMinusTwo ~ dgamma( 0.01 , 0.01 )  # mean=1 , sd=10 (generic vague)
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Run the chains:

parameters = c( "thetaC" , "thetaT" , "multT" ,
                "moC" , "kappaC" , "multTMo" , "multTSD" )
adaptSteps = 1000
burnInSteps = 1000
numSavedSteps=20000
thinSteps=20
nChains = 3

runJagsOut <- run.jags( method=c("rjags","parallel")[2] ,
                        model="TEMPmodel.txt" , 
                        monitor=parameters , 
                        data=dataList ,  
                        #inits=initsList , 
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps , 
                        sample=ceiling(numSavedSteps/nChains) ,
                        thin=thinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )
mcmcMat = as.matrix(codaSamples)

# Examine the chains:
# Convergence diagnostics:
diagMCMC( codaObject=codaSamples , parName="moC" )
diagMCMC( codaObject=codaSamples , parName="kappaC" )
diagMCMC( codaObject=codaSamples , parName="multTMo" )
diagMCMC( codaObject=codaSamples , parName="multTSD" )

# # Posterior descriptives:

cexMain=1.25
cexMax = 2.5
cexMin = 1.0
cexSlope = (cexMax-cexMin)/(max(theData$ControlTotal)-min(theData$ControlTotal))
cexInter = cexMax - cexSlope*max(theData$ControlTotal)

#------------------------------------------------------------------------------
openGraph(height=2.75,width=5)
#layout(matrix(1:2,ncol=2))
par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
#plotPost( mcmcMat[,"moC"] , xlab=bquote("moC") , xlim=c(0,0.25) ,
#          main="Prop. Deaths in Control, Overall" , cex.main=cexMain )
plotPost( mcmcMat[,"multTMo"] , xlab=bquote(mu[meta]) , xlim=c(0.1,1.6) ,
          main="Multiplier for Treatment, Overall" ,  cex.main=cexMain ,
          compVal=1.0 , 
          ROPE=c(0.95,1.05) )
graphName = paste0(fileNameRoot,"-Overall")
saveGraph( file=graphName , type="eps" )
saveGraph( file=graphName , type="png" )

#------------------------------------------------------------------------------
openGraph(height=2.75,width=8)
layout(matrix(1:2,ncol=2))
par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
plotPost( mcmcMat[,"kappaC"] , xlab=bquote("kappaC") ,
          main="Consist. of Prop. Deaths, Control" , cex.main=cexMain )
plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
graphName = paste0(fileNameRoot,"-AcrossStudyVariance")
saveGraph( file=graphName , type="eps" )
saveGraph( file=graphName , type="png" )

#------------------------------------------------------------------------------
for ( sIdx in 1:length(theData$StudyID) ) {
  openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  cexPt = cexInter + cexSlope * theData$ControlTotal[sIdx]
#   plotPost( mcmcMat[,paste0("thetaC[",sIdx,"]")] , 
#             xlab=paste0("thetaC[",sIdx,"]") , xlim=c(0,0.25) ,
#             HDItextPlace=0.8 ,
#             main=paste0("Prop. Deaths in Control, Study ",sIdx) , cex.main=cexMain )
#   points( ControlProp[sIdx] , 0 , 
#           pch=24 , col="black" , cex = cexPt , lwd=2 , bg="grey" )
  plotPost( mcmcMat[,paste0("multT[",sIdx,"]")] , 
            xlab=bquote(mu[.(sIdx)]) , xlim=c(0.1,1.6) ,
            main=paste0("Multiplier for Treatment, Study ",sIdx) ,  
            cex.main=cexMain ,
            compVal=1.0 , 
            ROPE=c(0.95,1.05) )
  points(  multT[sIdx] , 0 , 
           pch=24 , col="black" , cex = cexPt , lwd=2 , bg="grey" )
  graphName = paste0(fileNameRoot,"-Study",sIdx)
  saveGraph( file=graphName , type="eps" )
  saveGraph( file=graphName , type="png" )
}

#------------------------------------------------------------------------------
# For comparing with Gelman et al., who used log odds:
logOdds = -( log( mcmcMat[,"moC"] / (1-mcmcMat[,"moC"]) ) 
             - 
               log( mcmcMat[,"moC"]*mcmcMat[,"multTMo"] 
                    / (1-mcmcMat[,"moC"]*mcmcMat[,"multTMo"]) ) )
openGraph(height=3.0,width=4)
#layout(matrix(1:2,ncol=2))
par( mar=c(3.5,0.5,3.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( logOdds , xlab="logOdds" , #xlim=c(0,0.25) ,
          main="Log Odds Ratio, Overall" , cex.main=cexMain )
graphName = paste0(fileNameRoot,"-LogOddsRatio")
saveGraph( file=graphName , type="eps" )
saveGraph( file=graphName , type="png" )

#------------------------------------------------------------------------------

