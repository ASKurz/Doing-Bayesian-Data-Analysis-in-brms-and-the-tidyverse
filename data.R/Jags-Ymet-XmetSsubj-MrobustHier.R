# Jags-Ymet-XmetSsubj-MrobustHier.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , yName="y" , sName="s" ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault) { 

  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  # Convert sName to consecutive integers:
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    s = s ,
    Nsubj = max(s)  # should equal length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
    Ntotal <- length(y)
    xm <- mean(x)
    ym <- mean(y)
    xsd <- sd(x)
    ysd <- sd(y)
    for ( i in 1:length(y) ) {
      zx[i] <- ( x[i] - xm ) / xsd
      zy[i] <- ( y[i] - ym ) / ysd
    }
  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      zy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] , 1/zsigma^2 , nu )
    }
    for ( j in 1:Nsubj ) {
      zbeta0[j] ~ dnorm( zbeta0mu , 1/(zbeta0sigma)^2 )  
      zbeta1[j] ~ dnorm( zbeta1mu , 1/(zbeta1sigma)^2 )
    }
    # Priors vague on standardized scale:
    zbeta0mu ~ dnorm( 0 , 1/(10)^2 )
    zbeta1mu ~ dnorm( 0 , 1/(10)^2 )
    zsigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    nu ~ dexp(1/30.0)
    # Transform to original scale:
    for ( j in 1:Nsubj ) {
      beta1[j] <- zbeta1[j] * ysd / xsd  
      beta0[j] <- zbeta0[j] * ysd  + ym - zbeta1[j] * xm * ysd / xsd 
    }
    beta1mu <- zbeta1mu * ysd / xsd  
    beta0mu <- zbeta0mu * ysd  + ym - zbeta1mu * xm * ysd / xsd 
    sigma <- zsigma * ysd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta1" , "beta0mu" , "beta1mu" ,
                  "zbeta0" , "zbeta1" ,  "zbeta0mu" , "zbeta1mu" ,
                  "zsigma", "sigma", "nu" , "zbeta0sigma" , "zbeta1sigma" )
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 2000
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

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" , sName="s" ,
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
  s = factor(data[,sName])
  nSubj = length(unique(s)) # should be same as max(s)
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  beta0mu = mcmcMat[,"beta0mu"]
  beta1mu = mcmcMat[,"beta1mu"]
  sigma = mcmcMat[,"sigma"]
  zbeta0mu = mcmcMat[,"zbeta0mu"]
  zbeta1mu = mcmcMat[,"zbeta1mu"]
  zsigma = mcmcMat[,"zsigma"]
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
    pairs( cbind( beta0mu , beta1mu , sigma , log10nu )[plotIdx,] ,
           labels=c( expression(mu[beta*0]) , expression(mu[beta*1]) , 
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
  openGraph(width=8,height=8)
  layout( matrix( 1:9 , nrow=3, byrow=TRUE ) )
  par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( beta0mu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta0 , ROPE=ropeBeta0 ,
                       xlab=bquote(mu[beta*0]) , main=paste("Intercept, Group Level") )
  histInfo = plotPost( beta1mu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta1 , ROPE=ropeBeta1 ,
                       xlab=bquote(mu[beta*1]) , main=paste("Slope, Group Level") )
  plot( beta1mu[plotIdx] , beta0mu[plotIdx] , 
        xlab=bquote(mu[beta*1]) , ylab=bquote(mu[beta*0]) ,
        col="skyblue" , cex.lab = 1.75 )
  histInfo = plotPost( zbeta0mu , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValBeta0 , ROPE=ropeBeta0 ,
                       xlab=bquote(zmu[beta*0]) , main=paste("Intercept, Group Level") )
  histInfo = plotPost( zbeta1mu , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValBeta1 , ROPE=ropeBeta1 ,
                       xlab=bquote(zmu[beta*1]) , main=paste("Slope, Group Level") )
  plot( zbeta1mu[plotIdx] ,zbeta0mu[plotIdx] , 
        xlab=bquote(zmu[beta*1]) , ylab=bquote(zmu[beta*0]) ,
        col="skyblue" , cex.lab = 1.75 )
  histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValSigma , ROPE=ropeSigma ,
                       xlab=bquote(sigma) , main=paste("Scale, Subj Level") )
  histInfo = plotPost( log10nu , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=NULL , ROPE=NULL ,
                       xlab=bquote(log10(nu)) , main=paste("Normality, Subj Level") )
  plot( log10nu[plotIdx] , sigma[plotIdx] , 
        xlab=bquote(log10(nu)) ,ylab=bquote(sigma) , 
        col="skyblue" , cex.lab = 1.75 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostMarg",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  nPanels=25
  nPlots = ceiling(nSubj/nPanels)
  for ( plotIdx in 1:nPlots ) {
    openGraph()
    par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
    layout(matrix(1:nPanels,nrow=5,byrow=TRUE))
    xRang = max(x)-min(x)
    yRang = max(y)-min(y)
    xLimMult = 0.1
    yLimMult = 0.1
    xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
    yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
    #for ( sIdx in unique(ceiling(seq(1,nSubj,length=nPanels))) ) { 
    for ( sIdx in ((plotIdx-1)*nPanels+1):min(nSubj,(plotIdx-1)*nPanels+nPanels)) { 
      thisSrows = (as.numeric(s)==sIdx)
      plot( x[thisSrows] , y[thisSrows] , 
            cex=1.0 , lwd=1 , col="black" , xlim=xLim , ylim=yLim ,
            xlab=xName , ylab=yName , cex.lab=1.0 ,
            main=paste0("Unit: ",levels(s)[sIdx]) , 
            cex.main=1.0  )
      # Superimpose a smattering of believable regression lines:
      nPredCurves=30
      xComb = seq(xLim[1],xLim[2],length=301)
      for ( i in floor(seq(1,chainLength,length=nPredCurves)) ) {
        b0 = mcmcMat[i,paste0("beta0[",sIdx,"]")]
        b1 = mcmcMat[i,paste0("beta1[",sIdx,"]")]
        lines( xComb , b0+b1*xComb , col="skyblue" )
      }
      points( x[thisSrows] , y[thisSrows] , pch=19 )
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPredSubj",plotIdx), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  openGraph()
  par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
  # Plot data values:
  xRang = max(x)-min(x)
  yRang = max(y)-min(y)
  xLimMult = 0.2
  yLimMult = 0.2
  xLim= c( min(x)-xLimMult*xRang , max(x)+xLimMult*xRang )
  yLim= c( min(y)-yLimMult*yRang , max(y)+yLimMult*yRang )
  plot( x , y , pch="" , cex=1.0 , col="black" , 
        xlim=xLim , ylim=yLim ,
        xlab=xName , ylab=yName , cex.lab=1.0 ,
        main="All Units" , cex.main=1.0  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=70
  for ( i in floor(seq(1,chainLength,length=nPredCurves)) ) {
    abline( mcmcMat[i,"beta0mu"] , mcmcMat[i,"beta1mu"] , col="skyblue" )
  }
  for ( sIdx in 1:nSubj ) {
    thisSrows = (as.numeric(s)==sIdx)
    lines( x[thisSrows] , y[thisSrows] , type="o" , pch=19 ) #, pch=sIdx , col=sIdx )
  }
  #
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPredAll",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
  # Individual parameter estimates (intercept, slope)
  beta0hdi = apply( mcmcMat[, grep("^beta0\\[",colnames(mcmcMat)) ] , 2 , HDIofMCMC )
  beta1hdi = apply( mcmcMat[, grep("^beta1\\[",colnames(mcmcMat)) ] , 2 , HDIofMCMC )
  beta0range = range( beta0hdi )
  beta1range = range( beta1hdi )
  library("car") # install.packages("car")
  openGraph()
  par( mar=c(2,2,1,0)+.5 , mgp=c(1.5,0.5,0) )
  plot( 0 , 0 , type="n" , 
        xlim=beta0range , ylim=beta1range ,
        xlab="Intercept" , ylab="Slope" , 
        main="Individual Parameters (mean and 95%CI)" )
  for ( s in 1:nSubj ) {
    b0s = mcmcMat[, grep(paste0("^beta0\\[",s,"\\]"),colnames(mcmcMat)) ]
    b1s = mcmcMat[, grep(paste0("^beta1\\[",s,"\\]"),colnames(mcmcMat)) ]
    dataEllipse(b0s,b1s,levels=c(0.95),plot.points=FALSE,col="skyblue",lwd=1)
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostSubjPar",sep=""), type=saveType)
  }
}

#===============================================================================
