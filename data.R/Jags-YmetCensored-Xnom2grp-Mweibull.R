# Censored Weibull model.
# June 2015, John K. Kruschke.

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
source("DBDA2E-utilities.R") # Must be in R's current working directory.
fileNameRoot="Jags-YmetCensored-Xnom2grp-Mweibull" # For output file names.
graphFileType="eps" # "eps","png",etc.

# In R, the dweibull is parameterized differently than in JAGS.
# R:    p(y) = dweibull( y ,    shape=a , scale=b      ) 
# R: weibullMode = b * ((a-1)/a)^(1/a)
# R from JAGS: b = (1/lambda)^(1/a)
# JAGS:   y  ~ dweib(        nu=shape=a , lambda=1/b^a )
# JAGS: weibullMode = (1/lambda)^(1/nu) * ((nu-1)/nu)^(1/nu)
# JAGS from R: lambda = 1/b^a

#..............................................................................
# Generate data from Weibull and censor it:

# Specify Weibull parameters for R:
b=c(50,60)   # Two components for two groups.
nGrp=length(b) # Number of groups, here assumed to be 2.
a=rep(9,nGrp) # Weibull parameterization in JAGS crashes when a<1.2 or so...
N=c(111,105) # N of each group.
fileNameRoot = paste0(fileNameRoot,"-",a[1],b[1],b[2],sum(N))
# Random data:
set.seed(47405)
x = c( rep(1,N[1]) , rep(2,N[2]) ) # Group code
y = NULL
for ( gIdx in 1:nGrp ) {
  y = c(y, rweibull(n=N[gIdx],shape=a[gIdx],scale=b[gIdx]) ) # R parameterization
}
# Censor the data:
# Choose thresholds so they play nice with histogram breaks:
thresh = c( mean(y)-0.5*sd(y) , mean(y) , mean(y)+0.75*sd(y) )
histBreaks = seq( mean(y)-10*sd(y) , mean(y)+10*sd(y) , 0.25*sd(y) ) 
lowBreak = min( which( histBreaks > min(y) ) ) - 1
highBreak = max( which( histBreaks < max(y) ) ) + 1
histBreaks = histBreaks[ lowBreak : highBreak ]
yOrig = y
b0 = ( y <  thresh[1] )
b1 = ( y >= thresh[1] & y < thresh[2] ) 
b2 = ( y >= thresh[2] & y < thresh[3] ) 
b3 = ( y >= thresh[3] ) 
yBin = y
yBin[b0] = 0
yBin[b1] = 1
yBin[b2] = 2
yBin[b3] = 3
y[ b1 ] = NA
y[ b3 ] = NA
threshMat = matrix( rep(thresh,length(y)) , nrow=length(y) , byrow=TRUE )

#..............................................................................
# From here down:
# JAGS will use y, yBin, threshMat, and x
# Specialized plots will use histBreaks and thresh


# Plot data with generating distribution:
openGraph(width=6,height=7)
layout( matrix(1:2,nrow=2) )
par( mar=c(3.5,3.5,4.0,1.0) , mgp=c(2.0,0.7,0) )
nu=a         # Merely for translating to JAGS parameterization
lambda=1/b^a # Merely for translating to JAGS parameterization
#weibullMode = (1/lambda)^(1/nu) * ((nu-1)/nu)^(1/nu)
weibullMode = b * ((a-1)/a)^(1/a)
for ( gIdx in 1:nGrp ) {
  histInfo = hist( y[x==gIdx] , probability=TRUE , breaks=histBreaks ,
                   xlab="y (data)" , ylab="Probability density" ,
                   cex.lab=1.5 , cex.main=1.25 ,
                   main=bquote(list( "Group "*.(gIdx) ,
                                     "Weibull:" ,
                                     a==.(a[gIdx]) , 
                                     b==.(b[gIdx]) , 
                                     lambda==.(signif(lambda[gIdx],3)) ,
                                     mode==.(signif(weibullMode[gIdx],3)) )) ,
                   col="red" , border="white" )
  text( max(histInfo$breaks) , max(histInfo$density) ,
        bquote( atop( N[censored]==.(sum(is.na(y[x==gIdx]))) ,
                      N[total]==.(length(y[x==gIdx])) ) ) ,
        adj=c(1,1) , cex=1.5 )
  text( x=mean(c(thresh[1],thresh[2])) , y=0 , adj=c(0.5,0) ,
        labels=bquote(.(round(100*sum(yBin==1 & x==gIdx)/sum(x==gIdx),0))*"%") )
  text( x=thresh[3] , y=0 , adj=c(-0.2,0) ,
        labels=bquote(.(round(100*sum(yBin==3 & x==gIdx)/sum(x==gIdx),0))*"%") )
  yGrid = seq(min(histInfo$breaks),max(histInfo$breaks),length=501)
  lines( yGrid , dweibull(yGrid,a[gIdx],b[gIdx])
         *length(y[x==gIdx])/sum(is.finite(y[x==gIdx])) , lwd=3 )
  abline( v=weibullMode[gIdx] , lty="dashed" )
}
saveGraph( file=paste0(fileNameRoot,"-Data") , type=graphFileType )

#..............................................................................

# Package data for JAGS:
dataList = list(
  y = y ,
  x = x , # JAGS model assumes x has consecutive integers starting at 1
  nGrp = nGrp ,
  yBin = yBin ,
  threshMat = threshMat ,
  Ntotal = length(y) #,
  #meanY = mean(y,na.rm=TRUE) ,
  #sdY = sd(y,na.rm=TRUE)
)

# Define the JAGS model:
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dweib( nu , lambda[x[i]] )  # JAGS parameterization
    yBin[i] ~ dinterval( y[i] , threshMat[i, ] )
  }
  nu <- a   # a is shape in R dweibull. Here assumed same for all groups.
  a ~ dgamma( 1.105 , 0.105) # mode 1, sd 10
  for ( g in 1:nGrp ) {
    lambda[g] <- 1/b[g]^a  # b is scale in R dweibull. Different for each group.
    weibMode[g] <- (1/lambda[g])^(1/nu) * ((nu-1)/nu)^(1/nu)
    b[g] ~ dgamma( 1.105 , 0.105) # mode 1, sd 10
  }
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Initialize the missing y values:
# intial values of censored data:
yInit = rep( NA , length(y) )
for ( i in 1:length(y) ) {
  if ( is.na(y[i]) ) { # if y is censored 
    if ( yBin[i]==0 ) { 
      yInit[i] = threshMat[i,1]-1 
    } else if ( yBin[i]==ncol(threshMat) ) { 
      yInit[i] = threshMat[i,ncol(threshMat)]+1 
    } else {
      yInit[i] = (threshMat[i,yBin[i]]+threshMat[i,yBin[i]+1])/2
    }
  }
}
initsList = list( y=yInit )

# find censored y values to monitor:
yImp1 = paste0("y[",min( which( yBin==1 ) ),"]")
xImp1 = x[ min( which( yBin==1 ) ) ]
yImp3 = paste0("y[",min( which( yBin==3 ) ),"]")
xImp3 = x[ min( which( yBin==3 ) ) ]

# Run JAGS model:

nChains=3
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" , 
                        monitor=c("nu","lambda","weibMode","a","b",
                                  yImp1, yImp3) , 
                        data=dataList ,  
                        inits=initsList , 
                        n.chains=nChains ,
                        adapt=500 ,
                        burnin=500 , 
                        sample=ceiling(15000/nChains) ,
                        thin=5 ,
                        summarise=FALSE ,
                        plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
mcmcMat = as.matrix(codaSamples) # convert to matrix format

#..............................................................................
# Examine the chains:

# Convergence diagnostics:
for ( parName in c("a",
                   "b[1]","lambda[1]","weibMode[1]",
                   "b[2]","lambda[2]","weibMode[2]") ) {
  diagMCMC( codaObject=codaSamples , parName=parName , 
            saveName=paste0(fileNameRoot,"-") ,
            saveType=graphFileType )
}

# Posterior descriptives:
openGraph(height=7.5,width=10)
layout( matrix(1:12,nrow=3,byrow=FALSE) )
par( mar=c(3.5,3.75,4.0,1.0) , mgp=c(2.0,0.7,0) )
#
plotPost( codaSamples[,"a"] , main="a" , xlab="Param. Value" ,
          cex.main=1.5 , cex.lab=1.0 )
chainIdx = round(seq(1,nrow(mcmcMat),length=700))
plot( mcmcMat[chainIdx,"a"] , mcmcMat[chainIdx,"weibMode[1]"] , 
      xlab="a" , ylab="weibMode[1]" , cex.lab=2.0 , col="skyblue" )
plot( mcmcMat[chainIdx,"a"] , mcmcMat[chainIdx,"weibMode[2]"] , 
      xlab="a" , ylab="weibMode[2]" , cex.lab=2.0 , col="skyblue" )
#
plotPost( mcmcMat[,"b[2]"]/mcmcMat[,"b[1]"] , compVal=1.0 ,
          main="b[2] / b[1]" , xlab="Ratio" , cex.main=1.5 , cex.lab=1.0 )
plotPost( codaSamples[,"b[1]"] , 
          xlim=range(mcmcMat[,grep( "^b" , colnames(mcmcMat) , value=TRUE )]) ,
          main="b[1]" , xlab="Param. Value" , cex.main=1.5 , cex.lab=1.0 )
plotPost( codaSamples[,"b[2]"] , 
          xlim=range(mcmcMat[,grep( "^b" , colnames(mcmcMat) , value=TRUE )]) ,
          main="b[2]" , xlab="Param. Value" , cex.main=1.5 , cex.lab=1.0 )
#
plotPost( mcmcMat[,"lambda[2]"]/mcmcMat[,"lambda[1]"] , compVal=1.0 ,
          main="lambda[2] / lambda[1]" , xlab="Ratio" , cex.main=1.5 , cex.lab=1.0 )
plotPost( codaSamples[,"lambda[1]"] , 
          xlim=quantile( 
            mcmcMat[,grep( "^la" , colnames(mcmcMat) , value=TRUE )] ,
            probs=c(0,0.99) ) ,
          main="lambda[1]" , xlab="Param. Value" , cex.main=1.5 , cex.lab=1.0 )
plotPost( codaSamples[,"lambda[2]"] , 
          xlim=quantile( 
            mcmcMat[,grep( "^la" , colnames(mcmcMat) , value=TRUE )] ,
            probs=c(0,0.99) ) ,
          main="lambda[2]" , xlab="Param. Value" , cex.main=1.5 , cex.lab=1.0 )
#
plotPost( mcmcMat[,"weibMode[2]"]/mcmcMat[,"weibMode[1]"] , compVal=1.0 ,
          main="weibMode[2] / weibMode[1]" , xlab="Ratio" , cex.main=1.5 , cex.lab=1.0 )
plotPost( codaSamples[,"weibMode[1]"] , 
          xlim=range(mcmcMat[,grep( "^we" , colnames(mcmcMat) , value=TRUE )]) ,
          main="weibMode[1]" , xlab="Param. Value" , cex.main=1.5 , cex.lab=1.0 )
plotPost( codaSamples[,"weibMode[2]"] , 
          xlim=range(mcmcMat[,grep( "^we" , colnames(mcmcMat) , value=TRUE )]) ,
          main="weibMode[2]" , xlab="Param. Value" , cex.main=1.5 , cex.lab=1.0 )
#
saveGraph( file=paste0(fileNameRoot,"-Post") , type=graphFileType )

# Imputed y values:
openGraph(height=7,width=6)
layout( matrix(1:2,nrow=2) )
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples[,yImp1] , main=paste0("Censored yBin=1, Group ",xImp1) , 
          xlab="Imputed y value" , xlim=range(yOrig) ,
          cex.main=1.25, cex.lab=1.25 )
plotPost( codaSamples[,yImp3] , main=paste0("Censored yBin=3, Group ",xImp3) , 
          xlab="Imputed y value" , xlim=range(yOrig) ,
          cex.main=1.25, cex.lab=1.25 )
saveGraph( file=paste0(fileNameRoot,"-Post-ImputedY") , type=graphFileType )

# Plot data with posterior predictives:
# Plot data with generating distribution:
openGraph(width=6,height=7)
layout( matrix(1:2,nrow=2) )
par( mar=c(3.5,3.5,4.0,1.0) , mgp=c(2.0,0.7,0) )
for ( gIdx in 1:nGrp ) {
  histInfo = hist( y[x==gIdx] , probability=TRUE , breaks=histBreaks ,
                   xlab="y (data)" , ylab="Probability density" ,
                   cex.lab=1.5 , cex.main=1.75 ,
                   main=bquote(list( "Group "*.(gIdx) )) ,
                   col="red" , border="white" )
  text( max(histInfo$breaks) , max(histInfo$density) ,
        bquote( atop( N[censored]==.(sum(is.na(y[x==gIdx]))) ,
                      N[total]==.(length(y[x==gIdx])) ) ) ,
        adj=c(1,1) , cex=1.5 )
  text( x=mean(c(thresh[1],thresh[2])) , y=0 , adj=c(0.5,0) ,
        labels=bquote(.(round(100*sum(yBin==1 & x==gIdx)/sum(x==gIdx),0))*"%") )
  text( x=thresh[3] , y=0 , adj=c(-0.2,0) ,
        labels=bquote(.(round(100*sum(yBin==3 & x==gIdx)/sum(x==gIdx),0))*"%") )
  yGrid = seq(min(histInfo$breaks),max(histInfo$breaks),length=501)
  chainIdx = round(seq(1,nrow(mcmcMat),length=20))
  for ( stepIdx in chainIdx ) {
    lines( yGrid , 
           dweibull(yGrid,
                    mcmcMat[stepIdx,paste0("a")],
                    mcmcMat[stepIdx,paste0("b[",gIdx,"]")])
           *length(y)/sum(is.finite(y)) , 
           lwd=1 , col="skyblue" )
  }
}
saveGraph( file=paste0(fileNameRoot,"-Post-Pred") , type=graphFileType )

#..............................................................................
