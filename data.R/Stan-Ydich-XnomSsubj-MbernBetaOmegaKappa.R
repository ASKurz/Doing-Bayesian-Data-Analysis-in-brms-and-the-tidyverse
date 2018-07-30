# Stan-Ydich-XnomSsubj-MbernBetaOmegaKappa.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
source("DBDA2E-utilities.R")
#===============================================================================

genMCMC = function( data , sName="s" , yName="y" ,  
                    numSavedSteps=50000 , saveName=NULL , thinSteps=1 ) { 
  require(rjags)
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  # N.B.: This function expects the data to be a data frame, 
  # with one component named y being a vector of integer 0,1 values,
  # and one component named s being a factor of subject identifiers.
  y = as.numeric(data[,yName])
  s = as.numeric(data[,sName]) # ensures consecutive integer levels
  # Do some checking that data make sense:
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  Ntotal = length(y)
  Nsubj = length(unique(s))
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    s = s ,
    Ntotal = Ntotal ,
    Nsubj = Nsubj
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=1> Nsubj ;
    int<lower=1> Ntotal ;
    int<lower=0,upper=1> y[Ntotal] ;
    int<lower=1> s[Ntotal] ; // notice Ntotal not Nsubj
  }
  parameters {
    real<lower=0,upper=1> theta[Nsubj] ; // individual prob correct
    real<lower=0,upper=1> omega ;        // group mode
    real<lower=0> kappaMinusTwo ;        // group concentration minus two
  }
  transformed parameters {
    real<lower=0> kappa ;  
    kappa <- kappaMinusTwo + 2 ;
  }
  model {
    omega ~ beta( 1 , 1 ) ;
    kappaMinusTwo ~ gamma( 0.01 , 0.01 ) ; // mean=1 , sd=10 (generic vague)
    // kappaMinusTwo ~ gamma( 1.105125 , 0.1051249 ) ;  # mode=1 , sd=10 
    theta ~ beta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 ) ; // vectorized
    for ( i in 1:Ntotal ) {
      y[i] ~ bernoulli( theta[s[i]] ) ;
    }
  }
  " # close quote for modelString
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  initsList = function() {
    thetaInit = rep(0,Nsubj)
    for ( sIdx in 1:Nsubj ) { # for each subject
      includeRows = ( s == sIdx ) # identify rows of this subject
      yThisSubj = y[includeRows]  # extract data of this subject
      resampledY = sample( yThisSubj , replace=TRUE ) # resample
      thetaInit[sIdx] = sum(resampledY)/length(resampledY) 
    }
    thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
    meanThetaInit = mean( thetaInit )
    kappaInit = 100 # lazy, start high and let burn-in find better value
    return( list( theta=thetaInit , omega=meanThetaInit , 
                  kappaMinusTwo=kappaInit-2 ) )
  }
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta","omega","kappa") # The parameters to be monitored
  burnInSteps = 500            # Number of steps to burn-in the chains
  nChains = 4                  # nChains should be 2 or more for diagnostics 
  
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
                       thin = thinSteps ,
                       init = initsList ) # optional  
  # Or, accomplish above in one "stan" command; note stanDso is not separate.
  
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

smryMCMC = function(  codaSamples , compVal=0.5 , rope=NULL , 
                      diffIdVec=NULL , compValDiff=0.0 , ropeDiff=NULL , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  Ntheta = length(grep("theta",colnames(mcmcMat)))
  summaryInfo = NULL
  rowIdx = 0
  # overall omega:
  summaryInfo = rbind( summaryInfo , 
                       summarizePost( mcmcMat[,"omega"] ,
                                      compVal=compVal , ROPE=rope ) )
  rowIdx = rowIdx+1
  rownames(summaryInfo)[rowIdx] = "omega"
  # kappa:
  summaryInfo = rbind( summaryInfo , 
                       summarizePost( mcmcMat[,"kappa"] ,
                                      compVal=NULL , ROPE=NULL ) )
  rowIdx = rowIdx+1
  rownames(summaryInfo)[rowIdx] = "kappa"
  # individual theta's:
  for ( tIdx in 1:Ntheta ) {
    parName = paste0("theta[",tIdx,"]")
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # differences of individual theta's:
  if ( !is.null(diffIdVec) ) {
    Nidx = length(diffIdVec)
    for ( t1Idx in 1:(Nidx-1) ) {
      for ( t2Idx in (t1Idx+1):Nidx ) {
        parName1 = paste0("theta[",diffIdVec[t1Idx],"]")
        parName2 = paste0("theta[",diffIdVec[t2Idx],"]")
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

plotMCMC = function( codaSamples , data , compVal=0.5 , rope=NULL , 
                     diffIdVec=NULL , compValDiff=0.0 , ropeDiff=NULL , 
                     saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  # N.B.: This function expects the data to be a data frame, 
  # with one component named y being a vector of integer 0,1 values,
  # and one component named s being a factor of subject identifiers.
  y = data$y
  s = as.numeric(data$s) # converts character to consecutive integer levels
  # Now plot the posterior:
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  # Plot omega, kappa:
  openGraph(width=3.5*2,height=3.0)
  par( mfrow=c(1,2) )
  par( mar=c(3.5,3,1,1) , mgp=c(2.0,0.7,0) )
  postInfo = plotPost( mcmcMat[,"kappa"] , compVal=NULL , ROPE=NULL ,
                       xlab=bquote(kappa) , main="" , 
                       xlim=c( min(mcmcMat[,"kappa"]),
                               quantile(mcmcMat[,"kappa"],probs=c(0.990)) ) )
  postInfo = plotPost( mcmcMat[,"omega"] , compVal=compVal , ROPE=rope ,
                       xlab=bquote(omega) , main="Group Mode" ,
                       xlim=quantile(mcmcMat[,"omega"],probs=c(0.005,0.995)) )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostOmega",sep=""), type=saveType)
  }
  # Plot individual theta's and differences:
  if ( !is.null(diffIdVec) ) {
    Nidx = length(diffIdVec)
    openGraph(width=2.5*Nidx,height=2.0*Nidx)
    par( mfrow=c(Nidx,Nidx) )
    for ( t1Idx in 1:Nidx ) {
      for ( t2Idx in 1:Nidx ) {
        parName1 = paste0("theta[",diffIdVec[t1Idx],"]")
        parName2 = paste0("theta[",diffIdVec[t2Idx],"]")
        if ( t1Idx > t2Idx) {  
          # plot.new() # empty plot, advance to next
          par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
          nToPlot = 700
          ptIdx = round(seq(1,chainLength,length=nToPlot))
          plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , cex.lab=1.75 ,
                 xlab=parName2 , ylab=parName1 , col="skyblue" )
          abline(0,1,lty="dotted")
        } else if ( t1Idx == t2Idx ) {
          par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost( mcmcMat[,parName1] , cex.lab = 1.75 , 
                               compVal=compVal , ROPE=rope , cex.main=1.5 ,
                               xlab=parName1 , main="" )
          includeRows = ( s == diffIdVec[t1Idx] ) # rows of this subject in data
          dataPropor = sum(y[includeRows])/sum(includeRows) 
          points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
        } else if ( t1Idx < t2Idx ) {
          par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                               compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                               xlab=paste0(parName1,"-",parName2) , main="" , 
                               cex.lab=1.75 )
          includeRows1 = ( s == diffIdVec[t1Idx] ) # rows of this subject in data
          dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
          includeRows2 = ( s == diffIdVec[t2Idx] ) # rows of this subject in data
          dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
          points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
        }
      }
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostTheta",sep=""), type=saveType)
  }
}

#===============================================================================
