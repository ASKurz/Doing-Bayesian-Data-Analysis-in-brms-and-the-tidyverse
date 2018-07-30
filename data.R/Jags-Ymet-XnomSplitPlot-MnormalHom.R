
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm , yName , xBetweenName , xWithinName , xSubjectName ,
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName])
  xBetween = as.numeric(as.factor(datFrm[,xBetweenName]))
  xBetweenLvl =  levels(as.factor(datFrm[,xBetweenName]))
  NxBetweenLvl = length(unique(xBetween))
  xWithin = as.numeric(as.factor(datFrm[,xWithinName]))
  xWithinLvl =  levels(as.factor(datFrm[,xWithinName]))
  NxWithinLvl = length(unique(xWithin))
  xSubject = as.numeric(as.factor(datFrm[,xSubjectName]))
  xSubjectLvl =  levels(as.factor(datFrm[,xSubjectName]))
  NxSubjectLvl = length(unique(xSubject))
  # between-subject factor level of each subject:
  xBetweenOfSubject = aggregate( xBetween , list(xSubject) , FUN=unique )$x
  Ntotal = length(y)
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    xBetween = xBetween ,
    NxBetweenLvl = NxBetweenLvl ,
    xWithin = xWithin ,
    NxWithinLvl = NxWithinLvl ,
    xSubject = xSubject ,
    NxSubjectLvl = NxSubjectLvl ,
    xBetweenOfSubject = xBetweenOfSubject ,
    #xSubjectOfBetween = xSubjectOfBetween ,
    Ntotal = Ntotal ,
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD ,
    agammaShRa = agammaShRa 
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
    for ( r in 1:Ntotal ) {
      y[r] ~ dnorm( mu[r] , 1/sigma^2 )
      mu[r] <- ( a0 + aB[xBetween[r]] + aW[xWithin[r]] + aS[xSubject[r]] 
                 + aBxW[xBetween[r],xWithin[r]] )
    }
    sigma ~ dunif( ySD/100 , ySD*10 )
    a0 ~ dnorm( yMean , 1/(ySD*5)^2 )
    for ( i in 1:NxBetweenLvl ) { aB[i] ~ dnorm( 0.0 , 1/sigmaB^2 ) }
    sigmaB ~ dgamma(agammaShRa[1],agammaShRa[2])
    for ( j in 1:NxWithinLvl ) { aW[j] ~ dnorm( 0.0 , 1/sigmaW^2 ) }
    sigmaW ~ dgamma(agammaShRa[1],agammaShRa[2])
    for ( k in 1:NxSubjectLvl ) { aS[k] ~ dnorm( 0.0 , 1/sigmaS^2 ) }
    sigmaS ~ dgamma(agammaShRa[1],agammaShRa[2])
    for ( i in 1:NxBetweenLvl ) { for ( j in 1:NxWithinLvl ) {
      aBxW[i,j] ~ dnorm( 0.0 , 1/sigmaBxW^2 )
    } }
    sigmaBxW ~ dgamma(agammaShRa[1],agammaShRa[2])

    # For converting deflections to sum-to-zero:
    # JAGS doesn't like vectorized indices in mean(), so we use a 3D array.
    # Klunky, but it works.
    # Put predicted values into a 3D array, with empty cells set to zero:
    for ( k in 1:NxSubjectLvl ) { 
      for ( i in 1:NxBetweenLvl ) {
        for ( j in 1:NxWithinLvl ) {
          mSxBxW[k,i,j] <- ifelse( xBetweenOfSubject[k]==i ,
                                   a0 + aB[xBetweenOfSubject[k]] + aW[j] 
                                    + aBxW[xBetweenOfSubject[k],j] + aS[k] , 0 )
    } } }
    # Also put predicted values into a 2D array:
    for ( k in 1:NxSubjectLvl ) { for ( j in 1:NxWithinLvl ) {
      mSxW[k,j] <- ( a0 + aB[xBetweenOfSubject[k]] + aW[j] 
                     + aBxW[xBetweenOfSubject[k],j] + aS[k] )
    } }
    # Convert deflections to sum-to-zero:
    for ( i in 1:NxBetweenLvl ) { for ( j in 1:NxWithinLvl ) { 
      mBxW[i,j] <- ( sum( mSxBxW[1:NxSubjectLvl,i,j] ) 
                     / sum( mSxBxW[1:NxSubjectLvl,i,j]!=0 ) )
      bBxW[i,j] <- mBxW[i,j] - mB[i] - mW[j] + m
    } }
    for ( k in 1:NxSubjectLvl ) {
      mS[k] <- mean( mSxW[k,1:NxWithinLvl] )
      bS[k] <- mS[k] - mB[xBetweenOfSubject[k]]
    }
    m <- mean( mBxW[1:NxBetweenLvl,1:NxWithinLvl] )
    b0 <- m
    for ( i in 1:NxBetweenLvl ) {
      mB[i] <- mean( mBxW[i,1:NxWithinLvl] )
      bB[i] <- mB[i] - m
    }
    for ( j in 1:NxWithinLvl ) {
      mW[j] <- mean( mBxW[1:NxBetweenLvl,j] )
      bW[j] <- mW[j] - m
    }
        
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS it automatically...
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( # "a0" ,  "aB" ,  "aW" ,  "aBxW" , "aS" ,
    "mSxW" , "b0" ,  "bB" ,  "bW" ,  "bBxW" , "bS" , 
    "sigma" , "sigmaB" , "sigmaW" , "sigmaBxW" , "sigmaS"  )
  
  ## Using runjags:
  adaptSteps = 1000 
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
}

#===============================================================================

smryMCMC = function(  codaSamples , datFrm=NULL , 
                      xBetweenName=NULL , xWithinName=NULL , xSubjectName=NULL ,
                      xBetweenContrasts=NULL , 
                      xWithinContrasts=NULL , 
                      xBetweenWithinContrasts=NULL ,
                      saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(xBetweenName) & !is.null(xWithinName) 
       & !is.null(xSubjectName) ) {
    xBetweenLvl = levels(as.factor(datFrm[,xBetweenName]))
    xWithinLvl = levels(as.factor(datFrm[,xWithinName]))
    xSubjectLvl = levels(as.factor(datFrm[,xSubjectName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(xBetweenName) & !is.null(xWithinName) 
         & !is.null(xSubjectName) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "blah[12,34]" then pull out blah, 12 and 34:
      strparts = unlist( strsplit( parName , "\\[|,|\\]"  ) )
      # if there are only the param name and a single index:
      if ( length(strparts)==2 ) { 
        # if param name refers to Within:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="W" ) { 
          thisRowName = paste( thisRowName , xWithinLvl[as.numeric(strparts[2])] )
        }
        # if param name refers to factor Between:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="B" ) { 
          thisRowName = paste( thisRowName , xBetweenLvl[as.numeric(strparts[2])] )
        }
        # if param name refers to factor Subject:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="S" ) { 
          thisRowName = paste( thisRowName , xSubjectLvl[as.numeric(strparts[2])] )
        }
      }
      # if there are the param name and two indices:
      if ( length(strparts)==3 ) { 
        if ( substr(strparts[1],nchar(strparts[1])-2,nchar(strparts[1]))=="SxW" ) { 
          thisRowName = paste( thisRowName , 
                               xSubjectLvl[as.numeric(strparts[2])] , 
                               xWithinLvl[as.numeric(strparts[3])] )
        }
        if ( substr(strparts[1],nchar(strparts[1])-2,nchar(strparts[1]))=="BxW" ) { 
          thisRowName = paste( thisRowName , 
                               xBetweenLvl[as.numeric(strparts[2])] , 
                               xWithinLvl[as.numeric(strparts[3])] )
        }
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # Begin contrasts:
  if ( !is.null(xBetweenContrasts) | !is.null(xWithinContrasts) 
       | !is.null(xBetweenWithinContrasts) ) {
    if ( is.null(datFrm) | is.null(xBetweenName) | is.null(xWithinName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      
      xBetween = as.numeric(as.factor(datFrm[,xBetweenName]))
      xBetweenLvl = levels(as.factor(datFrm[,xBetweenName]))
      xWithin = as.numeric(as.factor(datFrm[,xWithinName]))
      xWithinLvl = levels(as.factor(datFrm[,xWithinName]))
      # xBetween contrasts:
      if ( !is.null(xBetweenContrasts) ) {
        for ( cIdx in 1:length(xBetweenContrasts) ) {
          thisContrast = xBetweenContrasts[[cIdx]]
          left = right = rep(FALSE,length(xBetweenLvl))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | xBetweenLvl==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | xBetweenLvl==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("bB[",1:length(xBetweenLvl),"]",sep="")] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste(thisContrast[[1]],collapse=""), ".v.",
                   paste(thisContrast[[2]],collapse=""),sep="") )
        }
      } # end between contrasts.
      # xWithin contrasts:
      if ( !is.null(xWithinContrasts) ) {
        for ( cIdx in 1:length(xWithinContrasts) ) {
          thisContrast = xWithinContrasts[[cIdx]]
          left = right = rep(FALSE,length(xWithinLvl))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | xWithinLvl==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | xWithinLvl==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("bW[",1:length(xWithinLvl),"]",sep="")] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste(thisContrast[[1]],collapse=""), ".v.",
                   paste(thisContrast[[2]],collapse=""),sep="") )
        }
      } # end within contrasts.
      # interaction contrasts:
      if ( !is.null(xBetweenWithinContrasts) ) {
        for ( cIdx in 1:length(xBetweenWithinContrasts) ) {
          thisContrast = xBetweenWithinContrasts[[cIdx]]
          # Between contrast:
          Bcontrast = thisContrast[[1]]
          BcontrastLeft = Bcontrast[[1]]
          BcontrastRight = Bcontrast[[2]]
          left = right = rep(FALSE,length(xBetweenLvl))
          for ( nIdx in 1:length( BcontrastLeft ) ) { 
            left = left | xBetweenLvl==BcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( BcontrastRight ) ) { 
            right = right | xBetweenLvl==BcontrastRight[nIdx]
          }
          right = normalize(right)
          BcontrastCoef = left-right
          # Within contrast:
          Wcontrast = thisContrast[[2]]
          WcontrastLeft = Wcontrast[[1]]
          WcontrastRight = Wcontrast[[2]]
          left = right = rep(FALSE,length(xWithinLvl))
          for ( nIdx in 1:length( WcontrastLeft ) ) { 
            left = left | xWithinLvl==WcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( WcontrastRight ) ) { 
            right = right | xWithinLvl==WcontrastRight[nIdx]
          }
          right = normalize(right)
          WcontrastCoef = left-right
          # interaction contrast:
          contrastCoef = cbind( as.vector( 
            outer( BcontrastCoef , WcontrastCoef ) ) )
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,5)=="bBxW["] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste( paste(thisContrast[[1]][[1]],collapse=""), 
                          ".v.",
                          paste(thisContrast[[1]][[2]],collapse=""),sep=""),
                   ".x.",
                   paste( paste(thisContrast[[2]][[1]],collapse=""), 
                          ".v.",
                          paste(thisContrast[[2]][[2]],collapse=""),sep=""),
                   sep="") )
        }
      } # end interaction contrasts
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
                     datFrm , yName="y" ,
                     xBetweenName=NULL , xWithinName=NULL , xSubjectName=NULL ,
                     xBetweenContrasts=NULL , 
                     xWithinContrasts=NULL , 
                     xBetweenWithinContrasts=NULL ,
                     saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  chainSub = round(seq(1,chainLength,length=20))
  y = datFrm[,yName]
  xBetween = as.numeric(as.factor(datFrm[,xBetweenName]))
  xBetweenLvl = levels(as.factor(datFrm[,xBetweenName]))
  xWithin = as.numeric(as.factor(datFrm[,xWithinName]))
  xWithinLvl = levels(as.factor(datFrm[,xWithinName]))
  xSubject = as.numeric(as.factor(datFrm[,xSubjectName]))
  xSubjectLvl = levels(as.factor(datFrm[,xSubjectName]))
  # Display data with posterior predictive distributions
  for ( Bidx in 1:length(xBetweenLvl) ) {
    openGraph(width=1.2*length(xWithinLvl),height=5)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    plot(-10,-10, 
         xlim=c(0.2,length(xWithinLvl)+0.1) , 
         xlab=paste(xWithinName,xBetweenName,sep="\n") , 
         xaxt="n" , ylab=yName ,
         ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) ,
         main="Data with Post. Pred.")
    axis( 1 , at=1:length(xWithinLvl) , tick=FALSE ,
          lab=paste( xWithinLvl , xBetweenLvl[Bidx] , sep="\n" ) )
    # plot subject data:
    SidxVec = unique( xSubject[ xBetween==Bidx ] )
    for ( Sidx in SidxVec ) {
      yVals = NULL
      for ( Widx in 1:length(xWithinLvl) ) {
        yVals = c( yVals , y[ xSubject==Sidx & xWithin==Widx ] )
      }
      plotCol="red"
      lines( 1:length(xWithinLvl)+runif(length(xWithinLvl),-0.05,0.05) , 
             yVals , type="o" , pch=1 , cex=1.5 , col=plotCol , lty="solid" )
      
    }
    # superimpose posterior predictive distrib:
    for ( Widx in 1:length(xWithinLvl) ) {
      for ( chnIdx in chainSub ) {
        m = ( mcmcMat[chnIdx,"b0"]
              + mcmcMat[chnIdx,paste0("bB[",Bidx,"]")]
              + mcmcMat[chnIdx,paste0("bW[",Widx,"]")]
              + mcmcMat[chnIdx,paste0("bBxW[",Bidx,",",Widx,"]")] )
        s = mcmcMat[chnIdx,"sigma"]
        normlim = qnorm( c(0.025,0.975) )
        yl = m+normlim[1]*s
        yh = m+normlim[2]*s
        ycomb=seq(yl,yh,length=201)
        ynorm = dnorm(ycomb,mean=m,sd=s)
        ynorm = 0.67*ynorm/max(ynorm)
        lines( Widx-ynorm , ycomb , col="skyblue" ) 
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",xBetweenLvl[Bidx]), type=saveType)
    }
    if ( length(xBetweenLvl) > 10 ) { dev.off() }
  } # end for Bidx
  
  #   # plot contrasts:
  if ( !is.null(xBetweenContrasts) | !is.null(xWithinContrasts) 
       | !is.null(xBetweenWithinContrasts) ) {
    if ( is.null(datFrm) | is.null(xBetweenName) | is.null(xWithinName) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      
      xBetween = as.numeric(as.factor(datFrm[,xBetweenName]))
      xBetweenLvl = levels(as.factor(datFrm[,xBetweenName]))
      xWithin = as.numeric(as.factor(datFrm[,xWithinName]))
      xWithinLvl = levels(as.factor(datFrm[,xWithinName]))
      
      # xBetween contrasts:
      if ( !is.null(xBetweenContrasts) ) {
        for ( cIdx in 1:length(xBetweenContrasts) ) {
          thisContrast = xBetweenContrasts[[cIdx]]
          left = right = rep(FALSE,length(xBetweenLvl))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | xBetweenLvl==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | xBetweenLvl==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("bB[",1:length(xBetweenLvl),"]",sep="")] 
                           %*% contrastCoef )
          
          openGraph(height=4,width=4)
          compInfo = plotPost( postContrast , xlab="Difference" ,
                               main=paste0( 
                                 paste(thisContrast[[1]],collapse="."), 
                                 "\nvs\n",
                                 paste(thisContrast[[2]],collapse=".") ),
                               compVal=thisContrast$compVal ,
                               ROPE=thisContrast$ROPE )
          saveGraph(file=paste0(fileNameRoot,
                                paste(thisContrast[[1]],collapse="."), 
                                "-vs-",
                                paste(thisContrast[[2]],collapse=".")),
                    type=graphFileType)
          
        }
      } # end between contrasts.
      
      # xWithin contrasts:
      if ( !is.null(xWithinContrasts) ) {
        for ( cIdx in 1:length(xWithinContrasts) ) {
          thisContrast = xWithinContrasts[[cIdx]]
          left = right = rep(FALSE,length(xWithinLvl))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | xWithinLvl==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | xWithinLvl==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("bW[",1:length(xWithinLvl),"]",sep="")] 
                           %*% contrastCoef )
          
          openGraph(height=4,width=4)
          compInfo = plotPost( postContrast , xlab="Difference" ,
                               main=paste0( 
                                 paste(thisContrast[[1]],collapse="."), 
                                 "\nvs\n",
                                 paste(thisContrast[[2]],collapse=".") ),
                               compVal=thisContrast$compVal ,
                               ROPE=thisContrast$ROPE )
          saveGraph(file=paste0(fileNameRoot,
                                paste(thisContrast[[1]],collapse="."), 
                                "-vs-",
                                paste(thisContrast[[2]],collapse=".")),
                    type=graphFileType)
          
        }
      } # end within contrasts.
      
      # interaction contrasts:
      if ( !is.null(xBetweenWithinContrasts) ) {
        for ( cIdx in 1:length(xBetweenWithinContrasts) ) {
          thisContrast = xBetweenWithinContrasts[[cIdx]]
          # Between contrast:
          Bcontrast = thisContrast[[1]]
          BcontrastLeft = Bcontrast[[1]]
          BcontrastRight = Bcontrast[[2]]
          left = right = rep(FALSE,length(xBetweenLvl))
          for ( nIdx in 1:length( BcontrastLeft ) ) { 
            left = left | xBetweenLvl==BcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( BcontrastRight ) ) { 
            right = right | xBetweenLvl==BcontrastRight[nIdx]
          }
          right = normalize(right)
          BcontrastCoef = left-right
          # Within contrast:
          Wcontrast = thisContrast[[2]]
          WcontrastLeft = Wcontrast[[1]]
          WcontrastRight = Wcontrast[[2]]
          left = right = rep(FALSE,length(xWithinLvl))
          for ( nIdx in 1:length( WcontrastLeft ) ) { 
            left = left | xWithinLvl==WcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( WcontrastRight ) ) { 
            right = right | xWithinLvl==WcontrastRight[nIdx]
          }
          right = normalize(right)
          WcontrastCoef = left-right
          # interaction contrast:
          contrastCoef = cbind( as.vector( 
            outer( BcontrastCoef , WcontrastCoef ) ) )
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,5)=="bBxW["] 
                           %*% contrastCoef )
          
          mainName = paste0(
            paste0(paste(thisContrast[[1]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[1]][[2]],collapse=".")),
            "\n(x)\n",
            paste0(paste(thisContrast[[2]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[2]][[2]],collapse=".")))
          openGraph(height=4,width=4)
          compInfo = plotPost( postContrast , xlab="Difference of Differences" ,
                               main=mainName ,
                               compVal=thisContrast$compVal ,
                               ROPE=thisContrast$ROPE )
          fileNameSuffix = paste0(
            paste0(paste(thisContrast[[1]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[1]][[2]],collapse=".")),
            "(x)",
            paste0(paste(thisContrast[[2]][[1]],collapse="."), 
                   ".v.",
                   paste(thisContrast[[2]][[2]],collapse=".")))
          saveGraph(file=paste0(fileNameRoot,fileNameSuffix),
                    type=graphFileType)
          
        }
      } # end interaction contrasts
      
      
    }
  }
  
}

# 
# #==============================================================================
