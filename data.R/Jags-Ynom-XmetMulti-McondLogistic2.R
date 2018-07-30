# Jags-Ynom-XmetMulti-McondLogistic.R 
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
  Nout = length( unique( y ) ) # should be same as max(y)
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  if ( Nout != 4 ) { stop("This model requires exactly 4 outcome categories.") }
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
    Nout = Nout
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
      y[i] ~ dcat( mu[1:Nout,i] ) 
      mu[1,i] <- phi[2,i] * phi[1,i]
      mu[2,i] <- (1-phi[2,i]) * phi[1,i]
      mu[3,i] <- phi[3,i] * (1-phi[1,i])
      mu[4,i] <- (1-phi[3,i]) * (1-phi[1,i])
      for ( r in 1:(Nout-1) ) {
        phi[r,i] <- ilogit( zbeta0[r] + sum( zbeta[r,1:Nx] * zx[i,1:Nx] ) )
      }
    }
    # Priors vague on standardized scale:
    for ( r in 1:(Nout-1) ) { 
      zbeta0[r] ~ dnorm( 0 , 1/20^2 )  
      for ( j in 1:Nx ) {
        zbeta[r,j] ~ dnorm( 0 , 1/20^2 )
      }
    }
    # Transform to original scale:
    for ( r in 1:(Nout-1) ) {
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
  openGraph(width=(Nx+1)*2.25,height=(Nout-1)*1.5)
  layout(matrix( 1:((Nx+1)*(Nout-1)) , byrow=TRUE , nrow=(Nout-1) ))
  par( mar=c(3.5,1,3,1) , mgp=c(2.0,0.7,0) )
  for ( oIdx in 1:(Nout-1) ) {
    # compute xlim range:
    xrange = NULL
    for ( r in 1:(Nout-1) ) {
      parName = paste0("beta0[",r,"]")
      xrange = range( c( xrange , mcmcMat[,parName] ) )
    }
    parName = paste0("beta0[",oIdx,"]")
    plotPost( mcmcMat[,parName] , xlab=parName , xlim=xrange ,
              border="skyblue" , HDItextPlace=0.9 ,
              main=paste0("Lambda: ",oIdx,", Pred: Intercept") )
    for ( xIdx in 1:Nx ) {
      # compute xlim range:
      xrange = NULL
      for ( r in 1:(Nout-1) ) {
        parName = paste0("beta[",r,",",xIdx,"]")
        xrange = range( c( xrange , mcmcMat[,parName] ) )
      }
      parName = paste0("beta[",oIdx,",",xIdx,"]")
      plotPost( mcmcMat[,parName] , xlab=parName , xlim=xrange , 
                border="skyblue" , HDItextPlace=0.9 ,
                main=paste0("Lambda: ",oIdx,", Pred: ",xIdx))
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
  # If exactly 2 predictors, then make X1,X2 plot:
  if ( Nx==2 ) {
    openGraph(width=6,height=7)
    par(mar=c(3,3,5.5,1),mgp=c(1.7,0.6,0),pty="s")
    # plot scatter of data:
    x1 = x[,1]
    x2 = x[,2]
    y = as.numeric(y)
    pCol = c("black","blue","red","green")
    plot( x1 , x2 , xlab=bquote(x[1]) , ylab=bquote(x[2]) , cex.lab=1.5 ,
          pch=as.character(y) , col=pCol[y] , cex=0.75 ,
          main=bquote( atop( 
            atop( {} ,
                  lambda["{"*1*","*2*"}|{"*1*","*2*","*3*","*4*"}"]  ) ,
            atop( lambda["{"*1*"}|{"*1*","*2*"}"] ,
                  lambda["{"*3*"}|{"*3*","*4*"}"]  ) ) ) , 
          cex.main=2.0 )
    # superimpose 50% logistic thresholds:
    chnVec = floor(seq(1,chainLength,length=20))
    for ( i in 1:(Nout-1) ) {
      for ( cIdx in chnVec ) {
        b0 = mcmcMat[cIdx,paste0("beta0[",i,"]")]
        b1 = mcmcMat[cIdx,paste0("beta[",i,",1]")]
        b2 = mcmcMat[cIdx,paste0("beta[",i,",2]")]
        if ( b2 != 0 ) {
          abline( -b0/b2 , -b1/b2 ,  col="skyblue" )
        } else {
          abline( v = -b0/b1 ,  col="skyblue" )
        }
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"Data",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
}
#===============================================================================
