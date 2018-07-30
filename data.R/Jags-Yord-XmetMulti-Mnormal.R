# Jags-Yord-XmetMulti-Mnormal.R 
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName , yName , 
                    numSavedSteps=10000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = as.matrix(data[,xName],ncol=length(xName))
  # Do some checking that data make sense:
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
  if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    warning("*** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS ***")
  }
  nYlevels = max(y)  
  thresh = rep(NA,nYlevels-1)
  thresh[1] = 1 + 0.5
  thresh[nYlevels-1] = nYlevels-1 + 0.5
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    nYlevels = nYlevels ,
    thresh = thresh ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1]
  )
  #
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
    for ( j in 1:Nx ) {
      xm[j]  <- mean(x[,j])
      xsd[j] <-   sd(x[,j])
      for ( i in 1:Ntotal ) {
        zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
      }
    }
  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dcat( pr[i,1:nYlevels] )
      pr[i,1] <- pnorm( thresh[1] , mu[i] , 1/sigma^2 )
      for ( k in 2:(nYlevels-1) ) {
        pr[i,k] <- max( 0 ,  pnorm( thresh[ k ] , mu[i] , 1/sigma^2 )
                           - pnorm( thresh[k-1] , mu[i] , 1/sigma^2 ) )
      }
      pr[i,nYlevels] <- 1 - pnorm( thresh[nYlevels-1] , mu[i] , 1/sigma^2 )
      mu[i] <- zbeta0 + sum( zbeta[1:Nx] * zx[i,1:Nx] ) 
    }
    # Priors vague on standardized scale:
    zbeta0 ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )
    for ( j in 1:Nx ) {
      zbeta[j] ~ dnorm( 0 , 1/(nYlevels)^2 )
    }
    zsigma ~ dunif( nYlevels/1000 , nYlevels*10 )
    # Transform to original scale:
    beta[1:Nx] <- ( zbeta[1:Nx] / xsd[1:Nx] )
    beta0 <- zbeta0  - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
    sigma <- zsigma
    for ( k in 2:(nYlevels-2) ) {  # 1 and nYlevels-1 are fixed
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
  parameters = c( "beta0" ,  "beta" ,  "sigma", "thresh" ,
                  "zbeta0" , "zbeta" , "zsigma" )
  adaptSteps = 500  # Number of steps to "tune" the samplers
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

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramName
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  zsigma = mcmcMat[,"zsigma"]
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  thresh  = mcmcMat[,grep("^thresh",colnames(mcmcMat))]
  sigma = mcmcMat[,"sigma"]
  #-----------------------------------------------------------------------------
#  # Compute R^2 for credible parameters:
#  YcorX = cor( y , x ) # correlation of y with each x predictor
#  Rsq = zbeta %*% matrix( YcorX , ncol=1 )
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
    pairs( cbind( beta0 , beta , sigma  )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , 
                     expression(sigma) ) , 
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
  openGraph(width=7.0,height=2.0)
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
  # If only a single predictor:
  if ( ncol(x)==1 ) {
    openGraph(width=7.0,height=5.0)
    par( mar=c(3.5,3.5,1,1) , mgp=c(2.25,0.7,0) )
    cexPt = 1.5
    plot( x , y , xlab=xName , ylab=yName , cex=cexPt , lwd=1.5 ,
          xlim=c( min(x)-0.2*(max(x)-min(x)) , max(x) ) ,
          ylim=c( min(y)-0.6 , max(y)+0.6 ) )
    # smattering of linear regression lines:
    nPredCurves = 20
    plotIdx = floor(seq(1,nrow(mcmcMat),length=nPredCurves))
    for ( i in plotIdx ) {
      abline( beta0[i] , beta[i,1] , col="skyblue" , lwd=2 )
    }
    # predictive categorical distributions:
    nSlice = 5
    curveXpos = seq(min(x),max(x),length=nSlice)
    curveWidth = (max(x)-min(x))/(nSlice+2)
    yComb = 1:max(y)
    for ( xIdx in 1:length(curveXpos) ) {   
      abline( v=curveXpos[xIdx] , col="grey" , lwd=2 ) # marks slice
      # rectangles of mean credible categorical distribution:
      pInt = pnorm ( colMeans(mcmcMat[,threshCols]) , 
                     mean = mean(beta0) + mean(beta)*curveXpos[xIdx] , 
                     sd = mean(sigma) )
      pInt = c(pInt,1) - c(0,pInt)
      xVals = curveWidth*pInt # /max(pInt)
      for ( yIdx in yComb ) {
        rect( xleft=curveXpos[xIdx]-xVals[yIdx] ,
              ybottom=yIdx-0.2 ,
              xright=curveXpos[xIdx] ,
              ytop=yIdx+0.4 , col="skyblue" , border="skyblue" )
      }
      # HDI segments:
      outProb=matrix(0,nrow=chainLength,ncol=max(y))
      for ( stepIdx in 1:chainLength ) {
        mu = mcmcMat[stepIdx,"beta0"] + mcmcMat[stepIdx,"beta"]*curveXpos[xIdx]
        threshCumProb = ( pnorm( ( mcmcMat[ stepIdx , 
                                            paste("thresh[",1:(max(y)-1),"]",sep="") ]
                                   - mu ) 
                                 / sigma[stepIdx] ) )
        outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
      }
      outHdi = apply( outProb , 2 , HDIofMCMC )
      outMean = apply( outProb , 2 , median , na.rm=TRUE )
      show(outMean)
      show(outHdi)
      for ( yIdx in 1:max(y) ) {
        segments( x0=curveXpos[xIdx]-curveWidth*outHdi[1,yIdx] , y0=yIdx+0.25 , 
                  x1=curveXpos[xIdx]-curveWidth*outHdi[2,yIdx] , y1=yIdx+0.25 , 
                  lwd=4 , col="darkgrey" )
      }
    }
    # re-plot data on top:
    points( x , y , cex=cexPt , lwd=1.5 )
    # save graph:
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"DataPostPred"), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # If exactly two predictors:
  if ( ncol(x)==2 ) {
    openGraph()
    par( mar=c(3.5,3.5,1,1) , mgp=c(2.25,0.7,0) )
    plot( x[,1] , x[,2] , pch=as.character(y) , xlab=xName[1] , ylab=xName[2] )
    nPredCurves = 3
    plotIdx = floor(seq(1,nrow(mcmcMat),length=nPredCurves))
    ltyCount=0
    for ( i in plotIdx ) {
      ltyCount = ltyCount+1
      for ( t in 1:ncol(thresh) ) {
        abline( (thresh[i,t]-beta0[i])/beta[i,2] , -beta[i,1]/beta[i,2] , 
                col="skyblue" , lwd=2 , lty=ltyCount )
      }
    }
    # save graph:
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"DataPostPred"), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=1 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(sigma) , main=paste("Std.Dev.") )
#   panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
#   histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
#                        xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMarg") )
  
  # Standardized scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
    histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zsigma , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*sigma) , main=paste("Std.Dev.") )
#   panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
#   histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
#                        xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}
#===============================================================================
