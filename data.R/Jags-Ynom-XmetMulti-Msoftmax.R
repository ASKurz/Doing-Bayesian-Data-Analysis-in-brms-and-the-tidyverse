# Jags-Ynom-XmetMulti-Msoftmax.R 
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
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = as.numeric(y) ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1] ,
    Nout = length( unique( y ) ) # should be same as max(y)
  )
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
      y[i] ~ dcat( explambda[1:Nout,i] ) # dcat normalizes its argument vector
      for ( r in 1:Nout ) {
        explambda[r,i] <- exp( zbeta0[r] + sum( zbeta[r,1:Nx] * zx[i,1:Nx] ) ) 
      }
    }
    # Priors vague on standardized scale:
    zbeta0[1] <- 0
    for ( j in 1:Nx ) { zbeta[1,j] <- 0 }
    for ( r in 2:Nout ) { # notice starts at 2
      zbeta0[r] ~ dnorm( 0 , 1/20^2 )  
      for ( j in 1:Nx ) {
        zbeta[r,j] ~ dnorm( 0 , 1/20^2 )
      }
    }
    # Transform to original scale:
    for ( r in 1:Nout ) {
      beta[r,1:Nx] <- zbeta[r,1:Nx] / xsd[1:Nx] 
      beta0[r] <- zbeta0[r] - sum( zbeta[r,1:Nx] * xm[1:Nx] / xsd[1:Nx] )
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
  parameters = c( "beta0" ,  "beta" ,  
                  "zbeta0" , "zbeta" )
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
  mcmcMat = as.matrix(codaSamples)
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

plotMCMC = function( codaSamples , data , xName , yName ,
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  Nx = ncol(x)
  Nout = length(unique(y))
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  openGraph(width=(Nx+1)*2.25,height=Nout*1.5)
  layout(matrix( 1:((Nx+1)*Nout) , byrow=TRUE , nrow=Nout ))
  par( mar=c(3.5,1,3,1) , mgp=c(2.0,0.7,0) )
  for ( oIdx in 1:Nout ) {
    # compute xlim range:
    xrange = NULL
    if ( oIdx != 1 ) {
      for ( r in 2:Nout ) {
        parName = paste0("beta0[",r,"]")
        xrange = range( c( xrange , mcmcMat[,parName] ) )
      }
    }
    parName = paste0("beta0[",oIdx,"]")
    plotPost( mcmcMat[,parName] , xlab=parName , xlim=xrange ,
              border="skyblue" , HDItextPlace=0.9 ,
              main=paste0("Out: ",oIdx,", Pred: Intercept") )
    for ( xIdx in 1:Nx ) {
      # compute xlim range:
      xrange = NULL
      if ( oIdx != 1 ) {
        for ( r in 2:Nout ) {
          parName = paste0("beta[",r,",",xIdx,"]")
          xrange = range( c( xrange , mcmcMat[,parName] ) )
        }
      }
      parName = paste0("beta[",oIdx,",",xIdx,"]")
      plotPost( mcmcMat[,parName] , xlab=parName , xlim=xrange , 
                border="skyblue" , HDItextPlace=0.9 ,
                main=paste0("Out: ",oIdx,", Pred: ",xIdx))
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
}
#===============================================================================
