# Stan-Ymet-Xmet-Mrobust.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , yName="y" , 
                    numSavedSteps=50000 , saveName=NULL ) { 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Ntotal = length(y) ,
    meanY = mean(y) ,
    sdY = sd(y) ,
    meanX = mean(x) ,
    sdX = sd(x)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=1> Ntotal ;
    real x[Ntotal] ;
    real y[Ntotal] ;
    real meanY ;
    real sdY ;
    real meanX ;
    real sdX ;
  }
  transformed data {
    real unifLo ;
    real unifHi ;
    real expLambda ;
    real beta0sigma ;
    real beta1sigma ;
    unifLo <- sdY/1000 ;
    unifHi <- sdY*1000 ;
    expLambda <- 1/30.0 ;
    beta1sigma <- 10*fabs(sdY/sdX) ;
    beta0sigma <- 10*fabs(meanX*sdY/sdX) ;
  }
  parameters {
    real beta0 ;
    real beta1 ;
    real<lower=0> nu ; 
    real<lower=0> sigma ; 
  }
  model {
    sigma ~ uniform( unifLo , unifHi ) ; 
    nu ~ exponential( expLambda ) ;
    beta0 ~ normal( 0 , beta0sigma ) ;
    beta1 ~ normal( 0 , beta1sigma ) ;
    for ( i in 1:Ntotal ) {
      y[i] ~ student_t( nu , beta0 + beta1 * x[i] , sigma  ) ;
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let Stan do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta1" ,  "sigma", "nu" )
  adaptSteps = 500  # Number of steps to "tune" the samplers
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

smryMCMC = function(  codaSamples , 
                      compValBeta0=NULL , ropeBeta0=NULL , 
                      compValBeta1=NULL , ropeBeta1=NULL , 
                      compValSigma=NULL , ropeSigma=NULL , 
                      saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "beta0" = summarizePost( mcmcMat[,"beta0"] , 
                                                compVal=compValBeta0 , 
                                                ROPE=ropeBeta0 ) )
  summaryInfo = rbind( summaryInfo , 
                       "beta1" = summarizePost( mcmcMat[,"beta1"] , 
                                                compVal=compValBeta1 , 
                                                ROPE=ropeBeta1 ) )
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
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     compValBeta0=NULL , ropeBeta0=NULL , 
                     compValBeta1=NULL , ropeBeta1=NULL , 
                     compValSigma=NULL , ropeSigma=NULL , 
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = data[,xName]
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  beta0 = mcmcMat[,"beta0"]
  beta1 = mcmcMat[,"beta1"]
  sigma = mcmcMat[,"sigma"]
  nu = mcmcMat[,"nu"]
  log10nu = log10(nu)
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
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
    pairs( cbind( beta0 , beta1 , sigma , log10nu )[plotIdx,] ,
           labels=c( expression(beta[0]) , expression(beta[1]) , 
                     expression(sigma) ,  expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  # Set up window and layout:
  nPtToPlot = 1000
  plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
  openGraph(width=8,height=5)
  layout( matrix( 1:6 , nrow=2, byrow=TRUE ) )
  par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta0 , ROPE=ropeBeta0 ,
                       xlab=bquote(beta[0]) , main=paste("Intercept") )
  histInfo = plotPost( beta1 , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta1 , ROPE=ropeBeta1 ,
                       xlab=bquote(beta[1]) , main=paste("Slope") )
  plot( beta1[plotIdx] , beta0[plotIdx] , 
        xlab=bquote(beta[1]) , ylab=bquote(beta[0]) ,
        col="skyblue" , cex.lab = 1.75 )
  histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValSigma , ROPE=ropeSigma ,
                       xlab=bquote(sigma) , main=paste("Scale") )
  histInfo = plotPost( log10nu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=NULL , ROPE=NULL ,
                       xlab=bquote(log10(nu)) , main=paste("Normality") )
  plot( log10nu[plotIdx] , sigma[plotIdx] , 
        xlab=bquote(log10(nu)) ,ylab=bquote(sigma) , 
        col="skyblue" , cex.lab = 1.75 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostMarg",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  openGraph()
  par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
  # Plot data values:
  postPredHDImass = 0.95
  xRang = max(x)-min(x)
  yRang = max(y)-min(y)
  xLimMult = 0.25
  yLimMult = 0.45
  xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
  yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
  plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
        xlab=xName , ylab=yName , cex.lab=1.5 ,
        main=paste( "Data w. Post. Pred. & ",postPredHDImass*100,"% HDI" ,sep="") , 
        cex.main=1.33  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=30
  xComb = seq(xLim[1],xLim[2],length=501)
  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    lines( xComb , beta0[i] + beta1[i]*xComb , col="skyblue" )
  }
  # Superimpose some vertical distributions to indicate spread:
  #source("HDIofICDF.R")
  nSlice = 5
  curveXpos = seq(min(x),max(x),length=nSlice)
  curveWidth = (max(x)-min(x))/(nSlice+2)
  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    for ( j in 1:length(curveXpos) ) {
      yHDI = HDIofICDF( qt , credMass=postPredHDImass , df=nu[i] )
      yComb = seq(yHDI[1],yHDI[2],length=75)
      xVals = dt( yComb , df=nu[i] ) 
      xVals = curveWidth * xVals / dt(0,df=nu[i])
      yPred = beta0[i] + beta1[i]*curveXpos[j] 
      yComb = yComb*sigma[i] + yPred
      lines( curveXpos[j] + xVals , yComb , col="skyblue" )
      lines( curveXpos[j] + 0*xVals , yComb , col="skyblue" , lwd=2 )
    }
  }
  # replot the data, in case they are obscured by lines:
  points( x , y , cex=1.5 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPred",sep=""), type=saveType)
  }
}

#===============================================================================
