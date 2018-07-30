# Jags-Yord-Xnom1grp-Mnormal.R 
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm, yName , numSavedSteps=50000 , thinSteps = 1,
                    saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = as.numeric(datFrm[,yName])
  # Do some checking that data make sense:
  if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
  if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    warning("*** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS ***")
  }
  Ntotal = length(y)
  # Threshold 1 and nYlevels-1 are fixed; other thresholds are estimated.
  # This allows all parameters to be interpretable on the response scale.
  nYlevels = max(y)  
  thresh = rep(NA,nYlevels-1)
  thresh[1] = 1 + 0.5
  thresh[nYlevels-1] = nYlevels-1 + 0.5
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    nYlevels = nYlevels ,
    thresh = thresh ,
    Ntotal = Ntotal 
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dcat( pr[i,1:nYlevels] )
      pr[i,1] <- pnorm( thresh[1] , mu , 1/sigma^2 )
      for ( k in 2:(nYlevels-1) ) {
        pr[i,k] <- max( 0 ,  pnorm( thresh[ k ] , mu , 1/sigma^2 )
                           - pnorm( thresh[k-1] , mu , 1/sigma^2 ) )
      }
      pr[i,nYlevels] <- 1 - pnorm( thresh[nYlevels-1] , mu , 1/sigma^2 )
    }
    mu ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )
    sigma ~ dunif( nYlevels/1000 , nYlevels*10 )
    for ( k in 2:(nYlevels-2) ) {  # 1 and nYlevels-1 are fixed, not stochastic
      thresh[k] ~ dnorm( k+0.5 , 1/2^2 )
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  #   # Initial values of MCMC chains based on data:
  #   muInit = c( mean(y[x==1]) , mean(y[x==2]) )
  #   sigmaInit = c( sd(y[x==1]) , sd(y[x==2]) )
  #   threshInit = 1:(nYlevels-1)+0.5
  #   threshInit[1] = NA
  #   threshInit[nYlevels-1] = NA
  #   # Regarding initial values in next line: (1) sigma will tend to be too big if 
  #   # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  #   # initial values keep the burn-in period moderate.
  #   initsList = list( mu=muInit, sigma=sigmaInit, nuMinusOne=4, thresh=threshInit )
  initsList = NULL
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma" , "thresh" )
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runJagsOut <- run.jags( method=runjagsMethod ,
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
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , compVal , #RopeEff=NULL , 
                      saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "mu" = summarizePost( mcmcMat[,"mu"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigma" = summarizePost( mcmcMat[,"sigma"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "effSz" = summarizePost( 
                         ( mcmcMat[,"mu"] - compVal ) / mcmcMat[,"sigma"] ,
                         compVal=compVal , ROPE=NULL ) )
  for ( colName in grep( "thresh" , colnames(mcmcMat) , value=TRUE ) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,colName] ) )
    rownames(summaryInfo)[nrow(summaryInfo)] = colName
  }
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , datFrm , yName , compVal , #RopeEff=NULL ,
                     showCurve=FALSE , pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  mu = mcmcMat[,"mu"]
  sigma = mcmcMat[,"sigma"]
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mu),by=length(mu)/nPtToPlot))
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
           labels=c( expression(mu) , expression(sigma) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Plot thresholds:
  threshCols = grep("thresh",colnames(mcmcMat),value=TRUE)
  threshMean = rowMeans( mcmcMat[,threshCols] )
  xLim = range(mcmcMat[,threshCols])
  nPtToPlot = 1000
  plotIdx = floor(seq(1,nrow(mcmcMat),length=nPtToPlot))
  openGraph(width=7.0,height=5.0)
  par( mar=c(3.5,3.5,1,1) , mgp=c(2.25,0.7,0) )
  plot( mcmcMat[plotIdx,threshCols[1]] , threshMean[plotIdx] , col="skyblue" ,
        xlim=xLim , xlab="Threshold" , ylab="Mean Threshold" )
  abline(v=mean(mcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="skyblue")
  for ( i in 2:length(threshCols) ) {
    points( mcmcMat[plotIdx,threshCols[i]] , threshMean[plotIdx] , col="skyblue" )
    abline(v=mean(mcmcMat[plotIdx,threshCols[i]]),lty="dashed",col="skyblue")
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Thresh"), type=saveType)
  }
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  openGraph(width=6.0,height=5.0)
  layout( matrix( c(1,2,3,4,5,6) , nrow=3, byrow=TRUE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  # Compute limits for plots of data with posterior pred. distributions
  y = as.numeric(datFrm[,yName])
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    warning("*** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS ***")
  }
  xLim = c( min(y)-0.5 , max(y)+0.5 )
  xBreaks = seq( xLim[1] , xLim[2] , 1 )  
  histInfo = hist(y,breaks=xBreaks,plot=FALSE)
  yMax = 1.2 * max( c( histInfo$density ) )
  xVec = seq( xLim[1] , xLim[2] , length=501 )
  #-----------------------------------------------------------------------------
  # 1. mu
  xlim = range( c( mu ) )
  histInfo = plotPost( mu , cex.lab = 1.75 ,
                       showCurve=showCurve , compVal=compVal ,
                       xlab=bquote(mu) , main=paste("Mean") , 
                       col="skyblue" )
  #-----------------------------------------------------------------------------
  # 2. Plot data y and smattering of posterior predictive curves:
  # Data histogram:
  histInfo = hist( y , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="pink" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data w. Post. Pred.") )
  # Posterior predictive probabilities of outcomes:
  outProb=matrix(0,nrow=chainLength,ncol=max(y))
  for ( stepIdx in 1:chainLength ) {
    threshCumProb = ( pnorm( ( mcmcMat[ stepIdx , 
                                        paste("thresh[",1:(max(y)-1),"]",sep="") ]
                               - mu[stepIdx] )
                             / sigma[stepIdx] ) )
    outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
  }
  outHdi = apply( outProb , 2 , HDIofMCMC )
  outMean = apply( outProb , 2 , median , na.rm=TRUE )
  show(outMean)
  points( x=1:max(y) , y=outMean  , pch=19 , cex=2 , col="skyblue" )
  segments( x0=1:max(y) , y0=outHdi[1,] , 
            x1=1:max(y) , y1=outHdi[2,] , lwd=4 , col="skyblue" )
  # annotate N:
  text( max(xVec) , yMax , bquote(N==.(length(y))) , adj=c(1.1,1.1) )
  #-----------------------------------------------------------------------------
  # 3. sigma
  xlim=range( c( sigma ) )
  histInfo = plotPost( sigma ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(sigma) , 
                       main=paste("Std. Dev.") , 
                       col="skyblue"  )
  #-----------------------------------------------------------------------------
  # 4. effect size. 
  effectSize = ( mu - compVal ) / sigma
  histInfo = plotPost( effectSize , compVal=0.0 , # ROPE=RopeEff ,
                       showCurve=showCurve ,
                       xlab=bquote( ( mu - .(compVal) ) / sigma ) ,
                       cex.lab=1.75 , main="Effect Size" ,
                       col="skyblue" )
  #-----------------------------------------------------------------------------  
  # 5. Thresholds:
  threshCols = grep("thresh",colnames(mcmcMat),value=TRUE)
  threshMean = rowMeans( mcmcMat[,threshCols] )
  xLim = range(mcmcMat[,threshCols])
  nPtToPlot = 500
  plotIdx = floor(seq(1,nrow(mcmcMat),length=nPtToPlot))
  plot( mcmcMat[plotIdx,threshCols[1]] , threshMean[plotIdx] , col="skyblue" ,
        xlim=xLim , xlab="Threshold" , ylab="Mean Threshold" )
  abline(v=mean(mcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="skyblue")
  for ( i in 2:length(threshCols) ) {
    points( mcmcMat[plotIdx,threshCols[i]] , threshMean[plotIdx] , col="skyblue" )
    abline(v=mean(mcmcMat[plotIdx,threshCols[i]]),lty="dashed",col="skyblue")
  }
  #   # Blank plot (was occupied by nu parameter in robust version)
  #   plot(1, ann=FALSE, axes=FALSE, xlim=c(0,1) , ylim=c(0,1) , 
  #        type="n" , xaxs="i" , yaxs="i" )  
  #   text(.5,.5,"[assumes y~normal]",adj=c(.5,.5))
  #-----------------------------------------------------------------------------
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
