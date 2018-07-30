# Jags-BivariateNormalScript.R 
# John Kruschke, November 2015.
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

# Load the functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.
want = c("ellipse")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

# Load the data:
myData = read.csv("HtWtData110.csv") # must have file in curr. work. dir.
# y must have named columns, with no missing values!
y = myData[,c("weight","height")] 

# Alternative data:
# myData = read.csv("HtWtData300.csv") 
# y = myData[,c("weight","height")] 
# y = y[1:30,] # optionally use subset, for watching width of HDI's 

# N = 100
# y1 = rnorm(N,10,2)
# y2 = 5+2*y1+rnorm(N,0,5)
# y3 = 50-2*y1+rnorm(N,0,5)
# y = cbind( Apple=y1 , Baker=y2 , Charlie=y3 )


# Standardize the data:
sdOrig = apply(y,2,sd)
meanOrig = apply(y,2,mean)
zy = apply(y,2,function(yVec){(yVec-mean(yVec))/sd(yVec)})
# Assemble data for sending to JAGS:
dataList = list(
  zy = zy ,
  Ntotal =  nrow(zy) ,
  Nvar = ncol(zy) ,
  # Include original data info for transforming to original scale:
  sdOrig = sdOrig ,
  meanOrig = meanOrig ,
  # For wishart (dwish) prior on inverse covariance matrix:
  zRscal = ncol(zy) ,  # for dwish prior
  zRmat = diag(x=1,nrow=ncol(zy))  # Rmat = diag(apply(y,2,var))
)

# Define the model:
modelString = "
model {
  for ( i in 1:Ntotal ) {
    zy[i,1:Nvar] ~ dmnorm( zMu[1:Nvar] , zInvCovMat[1:Nvar,1:Nvar] ) 
  }
  for ( varIdx in 1:Nvar ) { zMu[varIdx] ~ dnorm( 0 , 1/2^2 ) }
  zInvCovMat ~ dwish( zRmat[1:Nvar,1:Nvar] , zRscal )
  # Convert invCovMat to sd and correlation:
  zCovMat <- inverse( zInvCovMat )
  for ( varIdx in 1:Nvar ) { zSigma[varIdx] <- sqrt(zCovMat[varIdx,varIdx]) }
  for ( varIdx1 in 1:Nvar ) { for ( varIdx2 in 1:Nvar ) {
    zRho[varIdx1,varIdx2] <- ( zCovMat[varIdx1,varIdx2] 
                               / (zSigma[varIdx1]*zSigma[varIdx2]) )
  } }
  # Convert to original scale:
  for ( varIdx in 1:Nvar ) { 
    sigma[varIdx] <- zSigma[varIdx] * sdOrig[varIdx] 
    mu[varIdx] <- zMu[varIdx] * sdOrig[varIdx] + meanOrig[varIdx]
  }
  for ( varIdx1 in 1:Nvar ) { for ( varIdx2 in 1:Nvar ) {
    rho[varIdx1,varIdx2] <- zRho[varIdx1,varIdx2]
  } }
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )


# Run the chains:
require(rjags)
jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList ,  
                        n.chains=3 , n.adapt=500 )
update( jagsModel , n.iter=500 )
codaSamples = coda.samples( jagsModel , 
                            variable.names=c("mu","sigma","rho") ,
                            n.iter=10000 )

# Convergence diagnostics:
parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=codaSamples , parName=parName ) 
}


# Examine the posterior distribution:

mcmcMat = as.matrix(codaSamples)
chainLength = nrow(mcmcMat)
stepVec = floor(seq(1,chainLength,length=20)) # subset of chain
Nvar = ncol(y)

# Make plots of posterior distribution:

# Preparation -- define useful functions
library(ellipse)
expandRange = function( x , exMult=0.1 ) {
  lowVal = min(x)
  highVal = max(x)
  wid = max(x)-min(x)
  return( c( lowVal - exMult*wid , highVal + exMult*wid ) )
}

for ( varIdx in 1:Nvar ) {
  # Marginal posterior on means:
  openGraph()
  plotPost( mcmcMat[ , paste0("mu[",varIdx,"]") ] ,
            xlab=paste0("mu[",varIdx,"]") ,
            main=paste( "Mean of" , colnames(y)[varIdx] ) )
  # Marginal posterior on standard deviations:
  openGraph()
  plotPost( mcmcMat[ , paste0("sigma[",varIdx,"]") ] ,
            xlab=paste0("sigma[",varIdx,"]") ,
            main=paste( "SD of" , colnames(y)[varIdx] ) )
}

for ( varIdx1 in 1:(Nvar-1) ) {
  for ( varIdx2 in (varIdx1+1):Nvar ) {
    # Marginal posterior on correlation coefficient
    openGraph()
    plotPost( mcmcMat[ , paste0("rho[",varIdx1,",",varIdx2,"]") ] ,
              xlab=paste0("rho[",varIdx1,",",varIdx2,"]") ,
              main=paste( "Corr. of" , colnames(y)[varIdx1] , 
                          "and" , colnames(y)[varIdx2] ) )
    # Date with posterior ellipse
    openGraph()
    ellipseLevel = 0.90
    plot( y[,c(varIdx1,varIdx2)] , # pch=19 , 
          xlim=expandRange(y[,varIdx1]) , ylim=expandRange(y[,varIdx2]) ,
          xlab=colnames(y)[varIdx1] , ylab=colnames(y)[varIdx2] ,
          main=bquote("Data with posterior "*.(ellipseLevel)*" level contour") )
    # Posterior ellipses:
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx , 
                                paste0("rho[",varIdx1,",",varIdx2,"]") ] , 
                       scale=mcmcMat[ stepIdx , 
                                      c( paste0("sigma[",varIdx1,"]") ,
                                         paste0("sigma[",varIdx2,"]") ) ] , 
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[",varIdx1,"]") ,
                                          paste0("mu[",varIdx2,"]") ) ] , 
                       level=ellipseLevel ) , 
              type="l" , col="skyblue" , lwd=1 )
    }
    # replot data:
    points( y[,c(varIdx1,varIdx2)] )
  }
}


# Show data descriptives on console:
cor( y )
apply(y,2,mean)
apply(y,2,sd)

#---------------------------------------------------------------------------
