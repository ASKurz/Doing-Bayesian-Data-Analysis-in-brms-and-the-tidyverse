# Jags-YmetBinned-Xnom1grp-MnormalInterval.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data ,  yName , yBinName , threshVecName ,
                    numSavedSteps=50000 , saveName=NULL ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  ybin = data[,yBinName]
  threshMat = data[,threshVecName]
  if ( is.null(dim(threshMat)) ) { threshMat = cbind(threshMat) }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    ybin = ybin ,
    threshMat = threshMat ,
    Ntotal = Ntotal ,
    meanY = mean(y,na.rm=TRUE) ,
    sdY = sd(y,na.rm=TRUE)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm( mu , 1/sigma^2 )
      ybin[i] ~ dinterval( y[i] , threshMat[i, ] )
    }
    mu ~ dnorm( meanY , 1/(100*sdY)^2 )
    sigma ~ dunif( sdY/1000 , sdY*1000 )
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = mean(y,na.rm=TRUE) 
  sigma = sd(y,na.rm=TRUE) 
  # intial values of censored data:
  yInit = rep( NA , length(y) )
  for ( i in 1:length(y) ) {
    if ( is.na(y[i]) ) { # if y is censored 
      if ( ybin[i]==0 ) { 
        yInit[i] = threshMat[i,1]-1 
      } else if ( ybin[i]==ncol(threshMat) ) { 
        yInit[i] = threshMat[i,ncol(threshMat)]+1 
      } else {
        yInit[i] = (threshMat[i,ybin[i]]+threshMat[i,ybin[i]+1])/2
      }
    }
  }
  initsList = list( mu = mu , sigma = sigma , y=yInit )
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma")     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 4 
  thinSteps = 1
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , compValMu , # must specify compValMu
                      ropeMu=NULL , saveName=NULL ,
                      compValSigma=NULL , ropeSigma=NULL , 
                      compValEff=0.0 , ropeEff=c(-0.1,0.1) ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "mu" = summarizePost( mcmcMat[,"mu"] , 
                                             compVal=compValMu , ROPE=ropeMu ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigma" = summarizePost( mcmcMat[,"sigma"] , 
                                                compVal=compValSigma , 
                                                ROPE=ropeSigma ) )
  summaryInfo = rbind( summaryInfo , 
                       "effSz" = summarizePost( 
                         ( mcmcMat[,"mu"] - compValMu ) / mcmcMat[,"sigma"] ,
                         compVal=compValEff , ROPE=ropeEff ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , yName , yBinName , threshVecName , 
                     compValMu , # must specify compValMu
                     ropeMu=NULL , 
                     compValSigma=NULL , ropeSigma=NULL , 
                     compValEff=0.0 , ropeEff=NULL ,
                     showCurve=FALSE , pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  mu = mcmcMat[,"mu"]
  sigma = mcmcMat[,"sigma"]
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=5,height=5)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( mu , sigma )[plotIdx,] ,
           labels=c( expression(mu) , 
                     expression(sigma) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  openGraph(width=6.0,height=8.0*2.5/5)
  layout( matrix( c(2,3, 1,4) , nrow=2, byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  # Select thinned steps in chain for plotting of posterior predictive curves:
  nCurvesToPlot = 20
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  # Compute limits for plots of data with posterior pred. distributions
  y = data[,yName]
  threshMat = data[,threshVecName]
  if ( is.null(dim(threshMat)) ) { threshMat = cbind(threshMat) }
  xLim = c( min(y,na.rm=TRUE)-0.1*(max(y,na.rm=TRUE)-min(y,na.rm=TRUE)) , 
            max(y,na.rm=TRUE)+0.1*(max(y,na.rm=TRUE)-min(y,na.rm=TRUE)) )
  xBreaks = seq( xLim[1] , xLim[2] , 
                 length=ceiling((xLim[2]-xLim[1])/(sd(y,na.rm=TRUE)/4)) )
  histInfo = hist(y,breaks=xBreaks,plot=FALSE)
  yMax = 1.2 * max( histInfo$density )
  xVec = seq( xLim[1] , xLim[2] , length=501 )
  #-----------------------------------------------------------------------------
  # Plot data y and smattering of posterior predictive curves:
  histInfo = hist( y , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="red2" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , main="Data w. Post. Pred." )
  for ( t in 1:ncol(threshMat) ) {
    abline( v=threshMat[1,t] , lty="dashed" )
  }
  for ( stepIdx in 1:length(stepIdxVec) ) {
    lines(xVec, dnorm( xVec , mu[stepIdxVec[stepIdx]] , 
                       sigma[stepIdxVec[stepIdx]] ) * length(y)/sum(is.finite(y)) , 
          type="l" , col="skyblue" , lwd=1 )
  }
  text( max(xVec) , yMax , bquote(N==.(length(y))) , adj=c(1.1,1.1) )
  #-----------------------------------------------------------------------------
  histInfo = plotPost( mu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValMu , ROPE=ropeMu ,
                       xlab=bquote(mu) , main=paste("Mean") , 
                       col="skyblue" )
  #-----------------------------------------------------------------------------
  histInfo = plotPost( sigma , cex.lab=1.75 , showCurve=showCurve ,
                       compVal=compValSigma , ROPE=ropeSigma , 
                       xlab=bquote(sigma) , main=paste("Std. Dev.") , 
                       col="skyblue" )
  #-----------------------------------------------------------------------------
  effectSize = ( mu - compValMu ) / sigma
  histInfo = plotPost( effectSize , compVal=compValEff ,  ROPE=ropeEff ,
                       showCurve=showCurve , 
                       xlab=bquote( ( mu - .(compValMu) ) / sigma ),
                       cex.lab=1.75 , main="Effect Size" ,
                       col="skyblue" )
  #-----------------------------------------------------------------------------  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
