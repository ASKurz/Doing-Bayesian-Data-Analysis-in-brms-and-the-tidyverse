# Jags-Ymet-Xnom1grp-Mrobust.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , numSavedSteps=50000 , saveName=NULL ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    Ntotal = Ntotal ,
    meanY = mean(y) ,
    sdY = sd(y)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dt( mu , 1/sigma^2 , nu )
  }
  mu ~ dnorm( meanY , 1/(100*sdY)^2 )
  sigma ~ dunif( sdY/1000 , sdY*1000 )
  nu ~ dexp(1/30.0)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = mean(y) 
  sigma = sd(y) 
  initsList = list( mu = mu , sigma = sigma , nu = 5 )
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 4 
  thinSteps = 5
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
                       "nu" = summarizePost( mcmcMat[,"nu"] , 
                                             compVal=NULL , ROPE=NULL ) )
  summaryInfo = rbind( summaryInfo , 
                       "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) , 
                                             compVal=NULL , ROPE=NULL ) )
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

plotMCMC = function( codaSamples , data , compValMu , # must specify compValMu
                     ropeMu=NULL , 
                     compValSigma=NULL , ropeSigma=NULL , 
                     compValEff=0.0 , ropeEff=c(-0.1,0.1) ,
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
  nu = mcmcMat[,"nu"]
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7*3/5,height=7*3/5)
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
    pairs( cbind( mu , sigma , log10(nu) )[plotIdx,] ,
           labels=c( expression(mu) , 
                     expression(sigma) , 
                     expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  openGraph(width=6.0,height=8.0*3/5)
  layout( matrix( c(2,3,5, 1,4,6) , nrow=3, byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  # Select thinned steps in chain for plotting of posterior predictive curves:
  nCurvesToPlot = 20
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  # Compute limits for plots of data with posterior pred. distributions
  y = data
  xLim = c( min(y)-0.1*(max(y)-min(y)) , max(y)+0.1*(max(y)-min(y)) )
  xBreaks = seq( xLim[1] , xLim[2] , 
                 length=ceiling((xLim[2]-xLim[1])/(sd(y)/4)) )
  histInfo = hist(y,breaks=xBreaks,plot=FALSE)
  yMax = 1.2 * max( histInfo$density )
  xVec = seq( xLim[1] , xLim[2] , length=501 )
  #-----------------------------------------------------------------------------
  # Plot data y and smattering of posterior predictive curves:
  histInfo = hist( y , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="red2" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , main="Data w. Post. Pred." )
  for ( stepIdx in 1:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu[stepIdxVec[stepIdx]])/sigma[stepIdxVec[stepIdx]], 
                    df=nu[stepIdxVec[stepIdx]] )/sigma[stepIdxVec[stepIdx]] , 
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
                       compVal=compValSigma , ROPE=ropeSigma , cenTend="mode" ,
                       xlab=bquote(sigma) , main=paste("Scale") , 
                       col="skyblue" )
  #-----------------------------------------------------------------------------
  effectSize = ( mu - compValMu ) / sigma
  histInfo = plotPost( effectSize , compVal=compValEff ,  ROPE=ropeEff ,
                       showCurve=showCurve , cenTend="mode" ,
                       xlab=bquote( ( mu - .(compValMu) ) / sigma ),
                       cex.lab=1.75 , main="Effect Size" ,
                       col="skyblue" )
  #-----------------------------------------------------------------------------  
  postInfo = plotPost( log10(nu) , col="skyblue" , # breaks=30 ,
                       showCurve=showCurve ,
                       xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , 
                       cenTend="mode" ,
                       main="Normality" ) #  (<0.7 suggests kurtosis)
  #-----------------------------------------------------------------------------  
  # Blank plot 
  plot(1, ann=FALSE, axes=FALSE, xlim=c(0,1) , ylim=c(0,1) , 
       type="n" , xaxs="i" , yaxs="i" )  
  text(.5,.5,"[intentionally blank]",adj=c(.5,.5))
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
