
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , yName="y" , xNomName="xNom" , xMetName="xMet" ,
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic xNom,y variable names for model:
  y = as.numeric(datFrm[,yName])
  xNom = as.numeric(as.factor(datFrm[,xNomName]))
  xNomlevels = levels(as.factor(datFrm[,xNomName]))
  xMet = as.numeric(datFrm[,xMetName])
  Ntotal = length(y)
  NxNomLvl = length(unique(xNom))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  lmInfo = lm( datFrm[,yName] ~ datFrm[,xMetName] + datFrm[,xNomName] )
  residSD = sqrt(mean(lmInfo$residuals^2)) # residual root mean squared deviation
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    xNom = xNom ,
    xMet = xMet ,
    xMetMean = mean(xMet) ,
    Ntotal = Ntotal ,
    NxNomLvl = NxNomLvl ,
    # data properties for scaling the prior:
    xMetSD = sd(xMet) ,
    yMean = mean(y) ,
    ySD = sd(y) ,
    residSD = residSD ,
    agammaShRa = agammaShRa 
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm( mu[i] , 1/ySigma^2 )
      mu[i] <- a0 + a[xNom[i]] + aMet*( xMet[i] - xMetMean )
    }
    ySigma ~ dunif( residSD/100 , ySD*10 ) 
    a0 ~ dnorm( yMean , 1/(ySD*5)^2 )
    for ( j in 1:NxNomLvl ) { a[j] ~ dnorm( 0.0 , 1/aSigma^2 ) }
    aSigma ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
    aMet ~ dnorm( 0 , 1/(2*ySD/xMetSD)^2 ) 
    # Convert a0,a[] to sum-to-zero b0,b[] :
    b0 <- a0 + mean( a[1:NxNomLvl] ) + aMet*(-xMetMean)
    for ( j in 1:NxNomLvl ) { b[j] <- a[j] - mean( a[1:NxNomLvl] ) }
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  initsList = list(
    a0 = mean(y) - lmInfo$coefficients["datFrm[, xMetName]"] * mean(datFrm[,xMetName]) ,
    a = unname( c(
      lmInfo$coef["(Intercept)"] 
      + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] ,
      lmInfo$coef["(Intercept)"] 
      + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] 
      + lmInfo$coef[grep( "xNomName" , names(lmInfo$coef) )] ) ) - mean(y) ,
    aSigma = sd( c(
       lmInfo$coef["(Intercept)"] 
       + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] ,
       lmInfo$coef["(Intercept)"] 
       + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] 
       + lmInfo$coef[grep( "xNomName" , names(lmInfo$coef) )] ) ) ,
    ySigma = residSD ,
    aMet = unname(lmInfo$coefficients["datFrm[, xMetName]"]) 
    # Let JAGS do other parameters automatically...
  )
  #show( initsList ) 
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  
  parameters = c( "b0" ,  "b" , "aMet" , "aSigma" , "ySigma" )
  adaptSteps = 500 
  burnInSteps = 1000 
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  
#   nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
#   # Create, initialize, and adapt the model:
#   jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
#                           n.chains=nChains , n.adapt=adaptSteps )
#   # Burn-in:
#   cat( "Burning in the MCMC chain...\n" )
#   update( jagsModel , n.iter=burnInSteps )
#   # The saved MCMC chain:
#   cat( "Sampling final MCMC chain...\n" )
#   codaSamples = coda.samples( jagsModel , variable.names=parameters , 
#                               n.iter=nIter , thin=thinSteps )

  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}

#===============================================================================

smryMCMC = function(  codaSamples , datFrm=NULL , xNomName=NULL , xMetName=NULL ,
                      contrasts=NULL , saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(xNomName) ) {
    xNomlevels = levels(as.factor(datFrm[,xNomName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(xNomName) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "beta[12,34]" then pull out 12 and 34:
      levelVal = as.numeric( 
                  grep( "^[1-9]" , # grep only substrings that begin with digits.
                        # Return sll substrings split by "[" or "," or "]":
                        unlist( strsplit( parName , "\\[|,|\\]"  ) ) , 
                        value=TRUE ) )
      if ( length(levelVal) > 0 ) { 
        # Assumes there is only a single factor, i.e., levelVal has only entry: 
        thisRowName = paste(thisRowName,xNomlevels[levelVal]) 
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # All contrasts:
  if ( !is.null(contrasts) ) {
    if ( is.null(datFrm) | is.null(xNomName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      # contrasts:
      if ( !is.null(contrasts) ) {
        for ( cIdx in 1:length(contrasts) ) {
          thisContrast = contrasts[[cIdx]]
          left = right = rep(FALSE,length(xNomlevels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | xNomlevels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | xNomlevels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b[",1:length(xNomlevels),"]",sep="")] 
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
  }
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , datFrm , yName , xNomName , xMetName ,
                     contrasts=NULL , saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = as.numeric(datFrm[,yName])
  xNom = as.numeric(as.factor(datFrm[,xNomName]))
  xNomlevels = levels(as.factor(datFrm[,xNomName]))
  xMet = as.numeric(datFrm[,xMetName])
  Ntotal = length(y)
  NxNomLvl = length(unique(xNom))

  # Display data with posterior predictive distributions:
  for ( xNomLvlIdx in 1:NxNomLvl ) {
    # Open blank graph with appropriate limits:
    xLim = c( min(xMet)-0.2*(max(xMet)-min(xMet)) ,
              max(xMet)+0.2*(max(xMet)-min(xMet)) )
    yLim = c( min(y)-0.2*(max(y)-min(y)) , 
              max(y)+0.2*(max(y)-min(y)) )
    openGraph(width=4,height=5)
    par(mar=c(3,3,3,0.5)) # number of margin lines: bottom,left,top,right
    par(mgp=c(1.75,0.5,0)) # which margin lines to use for labels
    plot(2*max(xLim),2*max(yLim), # point out of range not seen 
         xlab=xMetName , xlim=xLim , ylab=yName , ylim=yLim , 
         main=paste(xNomlevels[xNomLvlIdx],"Data\nwith Post. Pred. Distrib.") ) 
    # plot credible regression lines and noise profiles:
    nSlice = 3
    curveXpos = seq(min(xMet),max(xMet),length=nSlice)
    curveWidth = (max(xMet)-min(xMet))/(nSlice+2)
    nPredCurves=30
    for ( i in floor(seq(from=1,to=nrow(mcmcMat),length=nPredCurves)) ) {
      intercept = mcmcMat[i,"b0"]+mcmcMat[i,paste0("b[",xNomLvlIdx,"]")]
      slope = mcmcMat[i,"aMet"]
      noise = mcmcMat[i,"ySigma"]
      abline( a=intercept , b=slope , col="skyblue" )
      for ( j in 1:nSlice ) {
        hdiLo = intercept+slope*curveXpos[j] - 1.96*noise
        hdiHi = intercept+slope*curveXpos[j] + 1.96*noise
        yComb = seq( hdiLo , hdiHi , length=75 )
        xVals = dnorm( yComb , mean=intercept+slope*curveXpos[j] , sd=noise ) 
        xVals = curveWidth * xVals / max(xVals)
        lines( curveXpos[j] - xVals , yComb , col="skyblue" )
        lines( curveXpos[j] - 0*xVals , yComb , col="skyblue" , lwd=2 )
      }
    }
    # plot data points:
    includeVec = ( xNom == xNomLvlIdx )
    xVals = xMet[includeVec]
    yVals = y[includeVec]
    points( xVals , yVals , pch=1 , cex=1.5 , col="red" )
    
    
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",xNomlevels[xNomLvlIdx]), 
                 type=saveType)
    }
  }
  
  # Display contrast posterior distributions:
  if ( !is.null(contrasts) ) {
    if ( is.null(datFrm) | is.null(xNomName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      for ( cIdx in 1:length(contrasts) ) {
        thisContrast = contrasts[[cIdx]]
        left = right = rep(FALSE,length(xNomlevels))
        for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
          left = left | xNomlevels==thisContrast[[1]][nIdx]
        }
        left = normalize(left)
        for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
          right = right | xNomlevels==thisContrast[[2]][nIdx]
        }
        right = normalize(right)
        contrastCoef = matrix( left-right , ncol=1 )
        postContrast = ( mcmcMat[,paste("b[",1:length(xNomlevels),"]",sep="")] 
                         %*% contrastCoef )
        openGraph(height=8,width=4)
        layout(matrix(1:2,ncol=1))
        plotPost( postContrast , xlab="Difference" ,
                  main=paste0( 
                    paste(thisContrast[[1]],collapse="."), 
                    "\nvs\n",
                    paste(thisContrast[[2]],collapse=".") ) ,
                  compVal=thisContrast$compVal , ROPE=thisContrast$ROPE )
        plotPost( postContrast/mcmcMat[,"ySigma"] , 
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

