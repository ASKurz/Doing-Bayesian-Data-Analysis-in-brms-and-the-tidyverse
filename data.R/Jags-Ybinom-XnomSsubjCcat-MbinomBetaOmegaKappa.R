# Jags-Ybinom-XnomSsubjCcat-MbinomBetaOmegaKappa.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
source("DBDA2E-utilities.R")
#===============================================================================

genMCMC = function( data , zName="z" , NName="N" , sName="s" , cName="c" ,
                    numSavedSteps=50000 , saveName=NULL , thinSteps=1 ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  require(rjags)
  require(runjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  # N.B.: This function expects the data to be a data frame, 
  # with one component z being a vector of integer # successes,
  # one component N being a vector of integer # attempts,
  # one component s being a factor of subject identifiers,
  # and one component c being a factor of category identifiers, with
  # subjects nested in categories.
  z = data[[zName]]
  N = data[[NName]]
  s = data[[sName]]
  c = data[[cName]]
  Nsubj = length(unique(s))
  Ncat =  length(unique(c))
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    z = z ,
    N = N ,
    c = as.numeric(c) , # c in JAGS is numeric, in R is possibly factor
    Nsubj = Nsubj ,
    Ncat = Ncat
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( sIdx in 1:Nsubj ) {
      z[sIdx] ~ dbin( theta[sIdx] , N[sIdx] )
      theta[sIdx] ~ dbeta( omega[c[sIdx]]*(kappa[c[sIdx]]-2)+1 , 
                           (1-omega[c[sIdx]])*(kappa[c[sIdx]]-2)+1 ) 
    }
    for ( cIdx in 1:Ncat ) {
      omega[cIdx] ~ dbeta( omegaO*(kappaO-2)+1 , 
                           (1-omegaO)*(kappaO-2)+1 )
      kappa[cIdx] <- kappaMinusTwo[cIdx] + 2
      kappaMinusTwo[cIdx] ~ dgamma( 0.01 , 0.01 ) # mean=1 , sd=10 (generic vague)
    }
    omegaO ~ dbeta( 1.0 , 1.0 ) 
    #omegaO ~ dbeta( 1.025 , 1.075 ) # mode=0.25 , concentration=2.1
    kappaO <- kappaMinusTwoO + 2
    kappaMinusTwoO ~ dgamma( 0.01 , 0.01 )  # mean=1 , sd=10 (generic vague)
    #kappaMinusTwoO ~ dgamma( 1.01005 , 0.01005012 )  # mode=1 , sd=100
    #kappaMinusTwoO ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10
    #kappaMinusTwoO ~ dgamma( 1.105125 , 0.01051249 )  # mode=10 , sd=100
  }
  " # close quote for modelString
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  initsList = function() {
    thetaInit = rep(NA,Nsubj)
    for ( sIdx in 1:Nsubj ) { # for each subject
      resampledZ = rbinom(1, size=N[sIdx] , prob=z[sIdx]/N[sIdx] )
      thetaInit[sIdx] = resampledZ/N[sIdx]
    }
    thetaInit = 0.001+0.998*thetaInit # keep away from 0,1    
    kappaInit = 100 # lazy, start high and let burn-in find better value
    return( list( theta=thetaInit , 
                  omega=aggregate(thetaInit,by=list(c),FUN=mean)$x ,
                  omegaO=mean(thetaInit) ,
                  kappaMinusTwo=rep(kappaInit-2,Ncat) ,
                  kappaMinusTwoO=kappaInit-2 ) )
  }
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta","omega","kappa","omegaO","kappaO") 
  adaptSteps = 500             # Number of steps to adapt the samplers
  burnInSteps = 500            # Number of steps to burn-in the chains
  
  useRunjags = TRUE
  if ( useRunjags ) {
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
  } else {
    # Create, initialize, and adapt the model:
    jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                            n.chains=nChains , n.adapt=adaptSteps )
    # Burn-in:
    cat( "Burning in the MCMC chain...\n" )
    update( jagsModel , n.iter=burnInSteps )
    # The saved MCMC chain:
    cat( "Sampling final MCMC chain...\n" )
    codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                                n.iter=ceiling(numSavedSteps*thinSteps/nChains), 
                                thin=thinSteps )
  }  
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , compVal=0.5 , rope=NULL , 
                      diffSVec=NULL , diffCVec=NULL , 
                      compValDiff=0.0 , ropeDiff=NULL , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = NULL
  rowIdx = 0
  # omega:
  for ( parName in grep("omega",colnames(mcmcMat),value=TRUE) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] ,
                                        compVal=compVal , ROPE=rope ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # kappa:
  for ( parName in grep("kappa",colnames(mcmcMat),value=TRUE) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] ,
                                        compVal=NULL , ROPE=NULL ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # theta:
  for ( parName in grep("theta",colnames(mcmcMat),value=TRUE) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] ,
                                        compVal=compVal , ROPE=rope ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # differences of theta's:
  if ( !is.null(diffSVec) ) {
    Nidx = length(diffSVec)
    for ( t1Idx in 1:(Nidx-1) ) {
      for ( t2Idx in (t1Idx+1):Nidx ) {
        parName1 = paste0("theta[",diffSVec[t1Idx],"]")
        parName2 = paste0("theta[",diffSVec[t2Idx],"]")
        summaryInfo = rbind( summaryInfo , 
                             summarizePost( mcmcMat[,parName1]-mcmcMat[,parName2] ,
                                            compVal=compValDiff , ROPE=ropeDiff ) )
        rowIdx = rowIdx+1
        rownames(summaryInfo)[rowIdx] = paste0(parName1,"-",parName2)
      }
    }
  }
  # differences of omega's:
  if ( !is.null(diffCVec) ) {
    Nidx = length(diffCVec)
    for ( t1Idx in 1:(Nidx-1) ) {
      for ( t2Idx in (t1Idx+1):Nidx ) {
        parName1 = paste0("omega[",diffCVec[t1Idx],"]")
        parName2 = paste0("omega[",diffCVec[t2Idx],"]")
        summaryInfo = rbind( summaryInfo , 
                             summarizePost( mcmcMat[,parName1]-mcmcMat[,parName2] ,
                                            compVal=compValDiff , ROPE=ropeDiff ) )
        rowIdx = rowIdx+1
        rownames(summaryInfo)[rowIdx] = paste0(parName1,"-",parName2)
      }
    }
  }
  # save:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  show( summaryInfo )
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , 
                     data , zName="z" , NName="N" , sName="s" , cName="c" ,
                     compVal=0.5 , rope=NULL , 
                     diffSList=NULL , diffCList=NULL , 
                     compValDiff=0.0 , ropeDiff=NULL , 
                     saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  # N.B.: This function expects the data to be a data frame, 
  # with one component z being a vector of integer # successes,
  # one component N being a vector of integer # attempts,
  # one component s being a factor of subject identifiers,
  # and one component c being a factor of category identifiers, with
  # subjects nested in categories.
  z = data[[zName]]
  N = data[[NName]]
  s = data[[sName]]
  c = data[[cName]]
  Nsubj = length(unique(s))
  Ncat =  length(unique(c))
  # Now plot the posterior:
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  # kappa:
  parNames = sort(grep("kappa",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 4
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,1,3.5,1) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.000,0.995))
  mainLab = c(levels(myData[[cName]]),"Overall")
  mainIdx = 0
  for ( parName in parNames ) {
    mainIdx = mainIdx+1
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=mainLab[mainIdx] , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Kappa",sep=""), type=saveType)
  }
  
  # omega:
  parNames = sort(grep("omega",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 4
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,1,3.5,1) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.001,0.999))
  mainLab = c(levels(myData[[cName]]),"Overall")
  mainIdx = 0
  for ( parName in parNames ) {
    mainIdx = mainIdx+1
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=mainLab[mainIdx] , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Omega",sep=""), type=saveType)
  }
  
  # Plot individual omega's and differences:
  if ( !is.null(diffCList) ) {
    for ( compIdx in 1:length(diffCList) ) {
      diffCVec = diffCList[[compIdx]]
      Nidx = length(diffCVec)
      temp=NULL
      for ( i in 1:Nidx ) {
        temp = c( temp , which(levels(myData[[cName]])==diffCVec[i]) )
      }
      diffCVec = temp
      openGraph(width=2.5*Nidx,height=2.0*Nidx)
      par( mfrow=c(Nidx,Nidx) )
      xLim = range(c( compVal, rope,
                      mcmcMat[,paste0("omega[",diffCVec,"]")] ))
      for ( t1Idx in 1:Nidx ) {
        for ( t2Idx in 1:Nidx ) {
          parName1 = paste0("omega[",diffCVec[t1Idx],"]")
          parName2 = paste0("omega[",diffCVec[t2Idx],"]")
          if ( t1Idx > t2Idx) {  
            # plot.new() # empty plot, advance to next
            par( mar=c(3,3,3,1) , mgp=c(2.0,0.7,0) , pty="s" )
            nToPlot = 700
            ptIdx = round(seq(1,chainLength,length=nToPlot))
            plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , 
                   cex.main=1.25 , cex.lab=1.25 , 
                   xlab=levels(myData[[cName]])[diffCVec[t2Idx]] , 
                   ylab=levels(myData[[cName]])[diffCVec[t1Idx]] , 
                   col="skyblue" )
            abline(0,1,lty="dotted")
          } else if ( t1Idx == t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1] , 
                                 compVal=compVal , ROPE=rope , 
                                 cex.main=1.25 , cex.lab=1.25 , 
                                 xlab=bquote(.(parName1)) ,
                                 main=levels(myData[[cName]])[diffCVec[t1Idx]] ,  
                                 xlim=xLim )
          } else if ( t1Idx < t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                                 compVal=compValDiff , ROPE=ropeDiff , 
                                 cex.main=1.25 , cex.lab=1.25 , 
                                 xlab=bquote("Difference of "*omega*"'s") , 
                                 main=paste( 
                                   levels(myData[[cName]])[diffCVec[t1Idx]] ,
                                   "-",
                                   levels(myData[[cName]])[diffCVec[t2Idx]] ) )
          }
        }
      }
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,"OmegaDiff",compIdx), type=saveType)
      }
    }
  }

  # Plot individual theta's and differences:
  if ( !is.null(diffSList) ) {
    for ( compIdx in 1:length(diffSList) ) {
      diffSVec = diffSList[[compIdx]]
      Nidx = length(diffSVec)
      temp=NULL
      for ( i in 1:Nidx ) {
        temp = c( temp , which(myData[[sName]]==diffSVec[i]) )
      }
      diffSVec = temp
      openGraph(width=2.5*Nidx,height=2.0*Nidx)
      par( mfrow=c(Nidx,Nidx) )
      xLim = range(c(compVal,rope,mcmcMat[,paste0("theta[",diffSVec,"]")],
                     z[diffSVec]/N[diffSVec]))
      for ( t1Idx in 1:Nidx ) {
        for ( t2Idx in 1:Nidx ) {
          parName1 = paste0("theta[",diffSVec[t1Idx],"]")
          parName2 = paste0("theta[",diffSVec[t2Idx],"]")
          if ( t1Idx > t2Idx) {  
            # plot.new() # empty plot, advance to next
            par( mar=c(3,3,3,1) , mgp=c(2.0,0.7,0) , pty="s" )
            nToPlot = 700
            ptIdx = round(seq(1,chainLength,length=nToPlot))
            plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , cex.lab=1.25 ,
                   xlab=s[diffSVec[t2Idx]] , 
                   ylab=s[diffSVec[t1Idx]] , 
                   col="skyblue" )
            abline(0,1,lty="dotted")
          } else if ( t1Idx == t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1] , 
                                 compVal=compVal , ROPE=rope , 
                                 cex.main=1.25 , cex.lab=1.25 , 
                                 xlab=bquote(.(parName1)) ,
                                 main=paste0( s[diffSVec[t1Idx]] , 
                                              " (",c[diffSVec[t1Idx]],")") ,  
                                 xlim=xLim )
            points( z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]] , 0 , 
                    pch="+" , col="red" , cex=3 )
            text( z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]] , 0 , 
                  bquote(list( z==.(z[diffSVec[t1Idx]]) ,
                               N==.(N[diffSVec[t1Idx]]) )) , 
                  adj=c( (z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]]-xLim[1])/
                           (xLim[2]-xLim[1]),-3.25) , col="red" )
          } else if ( t1Idx < t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                                 compVal=compValDiff , ROPE=ropeDiff , 
                                 cex.main=0.67 , cex.lab=1.25 , 
                                 xlab=bquote("Difference of "*theta*"'s") , 
                                 main=paste( 
                                   s[diffSVec[t1Idx]] , 
                                   " (",c[diffSVec[t1Idx]],")" ,
                                   "\n -",
                                   s[diffSVec[t2Idx]] , 
                                   " (",c[diffSVec[t2Idx]],")" ) )
            points( z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]]
                    - z[diffSVec[t2Idx]]/N[diffSVec[t2Idx]] , 0 , 
                    pch="+" , col="red" , cex=3 )
          }
        }
      }
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,"ThetaDiff",compIdx), type=saveType)
      }
    }
  }
}

#===============================================================================
