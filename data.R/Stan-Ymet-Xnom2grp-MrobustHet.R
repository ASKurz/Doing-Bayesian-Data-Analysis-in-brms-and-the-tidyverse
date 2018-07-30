# Stan-Ymet-Xnom2grp-MrobustHet.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm, yName="y" , xName="x" , numSavedSteps=50000 , 
                    saveName=NULL ) { 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))
  # Do some checking that data make sense:
  if ( any( x!=1 & x!=2 ) ) { stop("All x values must be 1 or 2.") }
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( length(x) != length(y) ) { stop("x and y must be same length.") }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    meanY = mean(y) ,
    sdY = sd(y)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=1> Ntotal ;
    int x[Ntotal] ;
    real y[Ntotal] ;
    real meanY ;
    real sdY ;
  }
  transformed data {
    real unifLo ;
    real unifHi ;
    real normalSigma ;
    real expLambda ;
    unifLo <- sdY/1000 ;
    unifHi <- sdY*1000 ;
    normalSigma <- sdY*100 ;
    expLambda <- 1/30.0 ;
  }
  parameters {
    real<lower=0> nu ; 
    real mu[2] ;               // 2 groups
    real<lower=0> sigma[2] ;   // 2 groups
  }
  model {
    sigma ~ uniform( unifLo , unifHi ) ; // vectorized
    mu ~ normal( meanY , normalSigma ) ; // vectorized
    nu ~ exponential( expLambda ) ;
    for ( i in 1:Ntotal ) {
      y[i] ~ student_t( nu , mu[x[i]] , sigma[x[i]] ) ;
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = c( mean(y[x==1]) , mean(y[x==2]) )
  sigma = c( sd(y[x==1]) , sd(y[x==2]) )
  # Regarding initial values in next line: (1) sigma will tend to be too big if 
  # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  # initial values keep the burn-in period moderate.
  initsList = list( mu = mu , sigma = sigma , nu = 5 )
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 4 
  thinSteps = 1
  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( model_code=modelString ) 
  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                       data = dataList , 
                       #pars = parameters , # optional
                       chains = nChains ,
                       iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                +burnInSteps ) , 
                       warmup = burnInSteps , 
                       #init = initsList , # optional
                       thin = thinSteps )
  # For consistency with JAGS-oriented functions in DBDA2E collection, 
  # convert stan format to coda format:
  codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                   function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
    save( stanFit , file=paste(saveName,"StanFit.Rdata",sep="") )
    save( stanDso , file=paste(saveName,"StanDso.Rdata",sep="") )
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
                       "nu" = summarizePost( mcmcMat[,"nu"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) ) )
  summaryInfo = rbind( summaryInfo , 
                       "effSz" = summarizePost( 
                         ( mcmcMat[,"mu[2]"]-mcmcMat[,"mu[1]"] ) 
                    / sqrt((mcmcMat[,"sigma[1]"]^2+mcmcMat[,"sigma[2]"]^2)/2) ,
                         compVal=0.0 , ROPE=RopeEff ) )
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
  nu = mcmcMat[,"nu"]
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
    pairs( cbind( mu1 , mu2 , sigma1 , sigma2 , log10(nu) )[plotIdx,] ,
           labels=c( expression(mu[1]) , expression(mu[2]) , 
                     expression(sigma[1]) , expression(sigma[2]) , 
                     expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  openGraph(width=6.0,height=8.0)
  layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  # Select thinned steps in chain for plotting of posterior predictive curves:
  nCurvesToPlot = 20
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  # Compute limits for plots of data with posterior pred. distributions
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))
  y1 = y[x==1]
  y2 = y[x==2]
  xLim = c( min(y)-0.1*(max(y)-min(y)) , max(y)+0.1*(max(y)-min(y)) )
  xBreaks = seq( xLim[1] , xLim[2] , 
                 length=ceiling((xLim[2]-xLim[1])/(mean(c(sd(y1),sd(y2)))/4)) )
  histInfo1 = hist(y1,breaks=xBreaks,plot=FALSE)
  histInfo2 = hist(y2,breaks=xBreaks,plot=FALSE)
  yMax = 1.2 * max( c( histInfo1$density , histInfo2$density ) )
  xVec = seq( xLim[1] , xLim[2] , length=501 )
  #-----------------------------------------------------------------------------
  # Plot data y1 and smattering of posterior predictive curves:
  histInfo = hist( y1 , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="red2" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data for",xLevels[1],"w. Post. Pred.") )
  for ( stepIdx in 1:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu1[stepIdxVec[stepIdx]])/sigma1[stepIdxVec[stepIdx]], 
                    df=nu[stepIdxVec[stepIdx]] )/sigma1[stepIdxVec[stepIdx]] , 
          type="l" , col="skyblue" , lwd=1 )
  }
  text( max(xVec) , yMax , bquote(N[1]==.(length(y1))) , adj=c(1.1,1.1) )
  #-----------------------------------------------------------------------------
  # Plot data y2 and smattering of posterior predictive curves:
  histInfo = hist( y2 , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="red2" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data for",xLevels[2],"w. Post. Pred.") )
  for ( stepIdx in 1:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu2[stepIdxVec[stepIdx]])/sigma2[stepIdxVec[stepIdx]], 
                    df=nu[stepIdxVec[stepIdx]] )/sigma2[stepIdxVec[stepIdx]] , 
          type="l" , col="skyblue" , lwd=1 )
  }
  text( max(xVec) , yMax , bquote(N[2]==.(length(y2))) , adj=c(1.1,1.1) )
  #-----------------------------------------------------------------------------  
  # Plot posterior distribution of parameter nu:
  histInfo = plotPost( log10(nu) , col="skyblue" , # breaks=30 ,
                       showCurve=showCurve ,
                       xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , 
                       cenTend="mode" ,
                       main="Normality" ) #  (<0.7 suggests kurtosis)
  
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
                       main=paste(xLevels[1],"Scale") , 
                       col="skyblue" , cenTend="mode" )
  histInfo = plotPost( sigma2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(sigma[2]) , 
                       main=paste(xLevels[2],"Scale") , 
                       col="skyblue" , cenTend="mode" )
  histInfo = plotPost( sigma2 - sigma1 , 
                       compVal=0 ,  showCurve=showCurve ,
                       xlab=bquote(sigma[2] - sigma[1]) , cex.lab = 1.75 , 
                       ROPE=RopeSdDiff ,
                       main="Difference of Scales" , col="skyblue" , 
                       cenTend="mode" )
  #-----------------------------------------------------------------------------
  # Plot effect size. 
  effectSize = ( mu2 - mu1 ) / sqrt( ( sigma1^2 + sigma2^2 ) / 2 )
  histInfo = plotPost( effectSize , compVal=0 ,  ROPE=RopeEff ,
                       showCurve=showCurve ,
                       xlab=bquote( (mu[2]-mu[1])
                                    /sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
                       cenTend="mode" , cex.lab=1.0 , main="Effect Size" ,
                       col="skyblue" )
  #-----------------------------------------------------------------------------
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
