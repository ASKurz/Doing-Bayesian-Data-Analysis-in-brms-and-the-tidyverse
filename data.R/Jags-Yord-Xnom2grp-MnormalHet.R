# Jags-Yord-Xnom2grp-MnormalHet.R 
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm, yName , xName , 
                    numSavedSteps=50000 , thinSteps = 1, saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))
  # Do some checking that data make sense:
  if ( any( x!=1 & x!=2 ) ) { stop("All x values must be 1 or 2.") }
  if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
  if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    cat("************************************************\n")
    cat("** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS **\n")
    cat("************************************************\n")
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
    x = x ,
    Ntotal = Ntotal 
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dcat( pr[i,1:nYlevels] )
      pr[i,1] <- pnorm( thresh[1] , mu[x[i]] , 1/sigma[x[i]]^2 )
      for ( k in 2:(nYlevels-1) ) {
        pr[i,k] <- max( 0 ,  pnorm( thresh[ k ] , mu[x[i]] , 1/sigma[x[i]]^2 )
                           - pnorm( thresh[k-1] , mu[x[i]] , 1/sigma[x[i]]^2 ) )
      }
      pr[i,nYlevels] <- 1 - pnorm( thresh[nYlevels-1] , mu[x[i]] , 1/sigma[x[i]]^2 )
    }
    for ( j in 1:2 ) { # 2 groups
      mu[j] ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )
      sigma[j] ~ dunif( nYlevels/1000 , nYlevels*10 )
    }
    for ( k in 2:(nYlevels-2) ) {  # 1 and nYlevels-1 are fixed, not stochastic
      thresh[k] ~ dnorm( k+0.5 , 1/2^2 )
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
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

smryMCMC = function(  codaSamples , RopeMuDiff=NULL , RopeSdDiff=NULL , 
                      RopeEff=NULL , saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "mu[1]" = summarizePost( mcmcMat[,"mu[1]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "mu[2]" = summarizePost( mcmcMat[,"mu[2]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "muDiff" = summarizePost( 
                         mcmcMat[,"mu[2]"] - mcmcMat[,"mu[1]"] , 
                         compVal=0.0 , ROPE=RopeMuDiff ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigma[1]" = summarizePost( mcmcMat[,"sigma[1]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigma[2]" = summarizePost( mcmcMat[,"sigma[2]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigmaDiff" = summarizePost( 
                         mcmcMat[,"sigma[2]"] - mcmcMat[,"sigma[1]"] , 
                         compVal=0.0 , ROPE=RopeSdDiff ) )
  summaryInfo = rbind( summaryInfo , 
                       "effSz" = summarizePost( 
                         ( mcmcMat[,"mu[2]"]-mcmcMat[,"mu[1]"] ) 
                    / sqrt((mcmcMat[,"sigma[1]"]^2+mcmcMat[,"sigma[2]"]^2)/2) ,
                         compVal=0.0 , ROPE=RopeEff ) )
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

plotMCMC = function( codaSamples , datFrm , yName="y" , xName="x" ,
                     RopeMuDiff=NULL , RopeSdDiff=NULL , RopeEff=NULL , 
                     showCurve=FALSE , pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # RopeMuDiff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of means.
  # RopeSdDiff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of standard deviations.
  # RopeEff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the effect size.
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  mu1 = mcmcMat[,"mu[1]"]
  mu2 = mcmcMat[,"mu[2]"]
  sigma1 = mcmcMat[,"sigma[1]"]
  sigma2 = mcmcMat[,"sigma[2]"]
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mu1),by=length(mu1)/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( mu1 , sigma1 , mu2 , sigma2 )[plotIdx,] ,
           labels=c( expression(mu[1]) , expression(sigma[1]) ,
                     expression(mu[2]) , expression(sigma[2]) ) , 
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
  openGraph(width=6.0,height=8.0)
  layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  # Compute limits for plots of data with posterior pred. distributions
  y = as.numeric(datFrm[,yName])
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    warning("*** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS ***")
  }
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))
  y1 = y[x==1]
  y2 = y[x==2]
  xLim = c( min(y)-0.5 , max(y)+0.5 )
  xBreaks = seq( xLim[1] , xLim[2] , 1 )  
  histInfo1 = hist(y1,breaks=xBreaks,plot=FALSE)
  histInfo2 = hist(y2,breaks=xBreaks,plot=FALSE)
  yMax = 1.2 * max( c( histInfo1$density , histInfo2$density ) )
  xVec = seq( xLim[1] , xLim[2] , length=501 )
  #-----------------------------------------------------------------------------
  # Plot data y1 and smattering of posterior predictive curves:
  # Data histogram:
  histInfo = hist( y1 , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="pink" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data for",xLevels[1],"w. Post. Pred.") )
  # Posterior predictive probabilities of outcomes:
  outProb=matrix(0,nrow=chainLength,ncol=max(y))
  for ( stepIdx in 1:chainLength ) {
    threshCumProb = ( pnorm( ( mcmcMat[ stepIdx , 
                                        paste("thresh[",1:(max(y)-1),"]",sep="") ]
                               - mu1[stepIdx] )
                             / sigma1[stepIdx] ) )
    outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
  }
  outHdi = apply( outProb , 2 , HDIofMCMC )
  outMean = apply( outProb , 2 , median , na.rm=TRUE )
  show(outMean)
  points( x=1:max(y) , y=outMean  , pch=19 , cex=2 , col="skyblue" )
  segments( x0=1:max(y) , y0=outHdi[1,] , 
            x1=1:max(y) , y1=outHdi[2,] , lwd=4 , col="skyblue" )
  # Annotate N:
  text( max(xVec) , yMax , bquote(N[1]==.(length(y1))) , adj=c(1.1,1.1) )
  #-----------------------------------------------------------------------------
  # Plot data y2 and smattering of posterior predictive curves:
  # Data histogram:
  histInfo = hist( y2 , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="pink" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data for",xLevels[2],"w. Post. Pred.") )
  # Posterior predictive probabilities of outcomes:
  outProb=matrix(0,nrow=chainLength,ncol=max(y))
  for ( stepIdx in 1:chainLength ) {
    threshCumProb = ( pnorm( ( mcmcMat[ stepIdx , 
                                        paste("thresh[",1:(max(y)-1),"]",sep="") ]
                               - mu2[stepIdx] )
                             / sigma2[stepIdx] ) )
    outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
  }
  outHdi = apply( outProb , 2 , HDIofMCMC )
  outMean = apply( outProb , 2 , median , na.rm=TRUE )
  show(outMean)
  points( x=1:max(y) , y=outMean  , pch=19 , cex=2 , col="skyblue" )
  segments( x0=1:max(y) , y0=outHdi[1,] , 
            x1=1:max(y) , y1=outHdi[2,] , lwd=4 , col="skyblue" )
  # Annotate N:
  text( max(xVec) , yMax , bquote(N[2]==.(length(y2))) , adj=c(1.1,1.1) )
  #-----------------------------------------------------------------------------  
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
  # Plot posterior distribution of parameters mu1, mu2, and their difference:
  xlim = range( c( mu1 , mu2 ) )
  histInfo = plotPost( mu1 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(mu[1]) , main=paste(xLevels[1],"Mean") , 
                       col="skyblue" )
  histInfo = plotPost( mu2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(mu[2]) , main=paste(xLevels[2],"Mean") , 
                       col="skyblue" )
  histInfo = plotPost( mu2-mu1 , compVal=0 ,  showCurve=showCurve ,
                       xlab=bquote(mu[2] - mu[1]) , cex.lab = 1.75 , 
                       ROPE=RopeMuDiff,
                       main="Difference of Means" , col="skyblue" )
  #-----------------------------------------------------------------------------
  # Plot posterior distribution of param's sigma1, sigma2, and their difference:
  xlim=range( c( sigma1 , sigma2 ) )
  histInfo = plotPost( sigma1 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(sigma[1]) , 
                       main=paste(xLevels[1],"Std. Dev.") , 
                       col="skyblue"  )
  histInfo = plotPost( sigma2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(sigma[2]) , 
                       main=paste(xLevels[2],"Std. Dev.") , 
                       col="skyblue"  )
  histInfo = plotPost( sigma2 - sigma1 , 
                       compVal=0 ,  showCurve=showCurve ,
                       xlab=bquote(sigma[2] - sigma[1]) , cex.lab = 1.75 , 
                       ROPE=RopeSdDiff ,
                       main="Difference of Std. Dev's" , col="skyblue"  )
  #-----------------------------------------------------------------------------
  # Plot effect size. 
  effectSize = ( mu2 - mu1 ) / sqrt( ( sigma1^2 + sigma2^2 ) / 2 )
  histInfo = plotPost( effectSize , compVal=0 ,  ROPE=RopeEff ,
                       showCurve=showCurve ,
                       xlab=bquote( (mu[2]-mu[1])
                                    /sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
                       cex.lab=1.0 , main="Effect Size" ,
                       col="skyblue" )
  #-----------------------------------------------------------------------------
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
