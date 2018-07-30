# Jags-Ymet-Xnom1fac-MrobustHet.R
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm , yName="y" , xName="x" ,
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xlevels = levels(as.factor(datFrm[,xName]))
  Ntotal = length(y)
  NxLvl = length(unique(x))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  # For hyper-prior on Factor 1 deflections:
  aeff = aggregate( y , list( x ) , mean )[,2] - yMean
  # For hyper-prior on deflections:
  aGammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # For lower limit on cell SDs:
  cellSDs = aggregate( y , list(x) , FUN=sd )
  cellSDs = cellSDs[ !is.na(cellSDs$x) ,]
  medianCellSD = median( cellSDs$x , na.rm=TRUE )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    NxLvl = NxLvl ,
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD ,
    aGammaShRa = aGammaShRa ,
    medianCellSD = medianCellSD
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dt(  a0 + a[x[i]] , 1/ySigma[x[i]]^2 , nu )
    }
    nu ~  dexp(1/30.0) 
    for ( j in 1:NxLvl ) { 
      ySigma[j] <- max( sigma[j] , medianCellSD/1000 ) # prevent zero ySigma
      sigma[j] ~ dgamma( ySigmaSh , ySigmaRa ) 
    }
    ySigmaSh <- 1 + ySigmaMode * ySigmaRa
    ySigmaRa <- ( ( ySigmaMode + sqrt( ySigmaMode^2 + 4*ySigmaSD^2 ) ) 
                  / ( 2*ySigmaSD^2 ) )
    ySigmaMode ~ dgamma( aGammaShRa[1] , aGammaShRa[2] ) 
    ySigmaSD ~ dgamma( aGammaShRa[1] , aGammaShRa[2] ) 
    #
    a0 ~ dnorm(yMean,1/(ySD*10)^2) 
    #
    for ( j in 1:NxLvl ) { a[j] ~ dnorm( 0.0 , 1/aSigma^2 ) }
    aSigma ~ dgamma( aGammaShRa[1] , aGammaShRa[2] ) 
  
    # Convert a0,a[] to sum-to-zero b0,b[] :
    for ( j in 1:NxLvl ) { m[j] <- a0 + a[j] } # cell means 
    b0 <- mean( m[1:NxLvl] )
    for ( j in 1:NxLvl ) { b[j] <- m[j] - b0 }
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it automatically...
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  # require(rjags)
  require(runjags)
  parameters = c( "b0" ,  "b" , "m" , "aSigma" ,
                  "ySigma" , "nu" , "ySigmaMode" , "ySigmaSD" )
  adaptSteps = 500 
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
}

#===============================================================================

smryMCMC = function(  codaSamples , datFrm=NULL , xName=NULL ,
                      contrasts=NULL , saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(xName) ) {
    xlevels = levels(as.factor(datFrm[,xName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(xName) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "beta[12,34]" then pull out 12 and 34:
      levelVal = as.numeric( 
        grep( "^[1-9]" , # grep only substrings that begin with digits.
              # Return sll substrings split by "[" or "," or "]":
              unlist( strsplit( parName , "\\[|,|\\]"  ) ) , 
              value=TRUE ) )
      if ( length(levelVal) > 0 ) { 
        # Assumes there is only a single factor, i.e., levelVal has only entry: 
        thisRowName = paste(thisRowName,xlevels[levelVal]) 
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # All contrasts:
  if ( !is.null(contrasts) ) {
    if ( is.null(datFrm) | is.null(xName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      for ( cIdx in 1:length(contrasts) ) {
        thisContrast = contrasts[[cIdx]]
        left = right = rep(FALSE,length(xlevels))
        for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
          left = left | xlevels==thisContrast[[1]][nIdx]
        }
        left = normalize(left)
        for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
          right = right | xlevels==thisContrast[[2]][nIdx]
        }
        right = normalize(right)
        contrastCoef = matrix( left-right , ncol=1 )
        postContrast = ( mcmcMat[,paste("b[",1:length(xlevels),"]",sep="")] 
                         %*% contrastCoef )
        summaryInfo = rbind( summaryInfo , 
                             summarizePost( postContrast ,
                                            compVal=thisContrast$compVal ,
                                            ROPE=thisContrast$ROPE ) )
        rownames(summaryInfo)[NROW(summaryInfo)] = (
          paste( paste(thisContrast[[1]],collapse=""), ".v.",
                 paste(thisContrast[[2]],collapse=""),sep="") )
      }
    }
  }
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================


plotMCMC = function( codaSamples , 
                     datFrm , yName="y" , xName="x" , contrasts=NULL ,
                     saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = datFrm[,yName]
  x = as.numeric(as.factor(datFrm[,xName]))
  xlevels = levels(as.factor(datFrm[,xName]))
  # Display data with posterior predictive distributions
  openGraph(width=min(10,1.25*length(xlevels)),height=5)
  plot(-1,0, 
       xlim=c(0.1,length(xlevels)+0.1) , 
       xlab=xName , xaxt="n" , ylab=yName ,
       ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) , 
       main="Data with Posterior Predictive Distrib.")
  axis( 1 , at=1:length(xlevels) , tick=FALSE , lab=xlevels )
  for ( xidx in 1:length(xlevels) ) {
    xPlotVal = xidx 
    yVals = y[ x==xidx ]
    points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) , 
            yVals , pch=1 , cex=1.5 , col="red" )
    chainSub = round(seq(1,chainLength,length=20))
    for ( chnIdx in chainSub ) {
      m = mcmcMat[chnIdx,paste("m[",xidx,"]",sep="")]
      s = mcmcMat[chnIdx,paste("ySigma[",xidx,"]",sep="")]
      nu = mcmcMat[chnIdx,"nu"]
      tlim = qt( c(0.025,0.975) , df=nu )
      yl = m+tlim[1]*s
      yh = m+tlim[2]*s
      ycomb=seq(yl,yh,length=201)
      #ynorm = dnorm(ycomb,mean=m,sd=s)
      #ynorm = 0.67*ynorm/max(ynorm)
      yt = dt( (ycomb-m)/s , df=nu )
      yt = 0.67*yt/max(yt)
      lines( xPlotVal-yt , ycomb , col="skyblue" ) 
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPred",sep=""), type=saveType)
  }
  if ( !is.null(contrasts) ) {
    if ( is.null(datFrm) | is.null(xName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      for ( cIdx in 1:length(contrasts) ) {
        thisContrast = contrasts[[cIdx]]
        left = right = rep(FALSE,length(xlevels))
        for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
          left = left | xlevels==thisContrast[[1]][nIdx]
        }
        left = normalize(left)
        for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
          right = right | xlevels==thisContrast[[2]][nIdx]
        }
        right = normalize(right)
        contrastCoef = matrix( left-right , ncol=1 )
        postContrast = ( mcmcMat[,paste("b[",1:length(xlevels),"]",sep="")] 
                         %*% contrastCoef )
        openGraph(height=8,width=4)
        layout(matrix(1:2,ncol=1))
        plotPost( postContrast , xlab="Difference" ,
                  main=paste0( 
                    paste(thisContrast[[1]],collapse="."), 
                    "\nvs\n",
                    paste(thisContrast[[2]],collapse=".") ) ,
                  compVal=thisContrast$compVal , ROPE=thisContrast$ROPE )
        plotPost( ( postContrast / sqrt( 
          ( mcmcMat[,paste0("ySigma[",1:length(xlevels),"]")]^2 
            %*% ( contrastCoef != 0 ) ) / sum(contrastCoef != 0 )) ) , 
          xlab="Effect Size" ,
          main=paste0( 
            paste(thisContrast[[1]],collapse="."), 
            "\nvs\n",
            paste(thisContrast[[2]],collapse=".") ) ,
          compVal=0.0 , 
          ROPE=c(-0.1,0.1) )
        
        if ( !is.null(saveName) ) {
          saveGraph( file=paste0(saveName, paste0( 
            paste(thisContrast[[1]],collapse=""), 
            ".v.",
            paste(thisContrast[[2]],collapse="") ) ), 
            type=saveType )
        }
      }
    }
  } # end if ( !is.null(contrasts) )
}
