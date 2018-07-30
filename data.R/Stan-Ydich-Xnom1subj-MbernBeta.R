# Stan-Ydich-Xnom1subj-MbernBeta.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
source("DBDA2E-utilities.R")
#===============================================================================

genMCMC = function( data , numSavedSteps=50000 , saveName=NULL ) { 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  if ( class(data)=="data.frame" ) {  # If data is a data.frame
    y = myData$y                      # then pull out the column named y
  } else {                            # else
    y = data                          # rename the data as y.
  }
  # Do some checking that data make sense:
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    Ntotal = Ntotal 
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=0> Ntotal ;
    int<lower=0,upper=1> y[Ntotal] ;
  }
  parameters {
    real<lower=0,upper=1> theta ;
  }
  model {
    theta ~ beta(1,1) ;
    y ~ bernoulli(theta) ; // implicitly vectorized
  }
  " # close quote for modelString
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  # Option 1: Use single initial value for all chains:
  #  thetaInit = sum(y)/length(y)
  #  initsList = list( theta=thetaInit )
  # Option 2: Use function that generates random values for each chain:
  initsList = function() {
    resampledY = sample( y , replace=TRUE )
    thetaInit = sum(resampledY)/length(resampledY)
    thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
    return( list( theta=thetaInit ) )
  }
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta")     # The parameters to be monitored
  burnInSteps = 500            # Stan defaults to iter/2 for overdispersed inits
  nChains = 4                  # nChains should be 2 or more for diagnostics 
  thinSteps = 4                # In Stan there is autocorrelation, so thin

#  stanCpp <- stanc( model_code = modelString ) # Translate to C++
#  stanDso <- stan_model( stanc_ret = stanCpp ) # Compile Stan DSO
  
  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( model_code=modelString ) 
  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                       data = dataList , 
                       pars = parameters , # optional
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

smryMCMC = function( codaSamples , compVal=NULL , rope=NULL , saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "theta" = summarizePost( mcmcMat[,"theta"] , 
                                             compVal=compVal , ROPE=rope ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  show( summaryInfo )
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , compVal=NULL , rope=NULL , 
                     saveName=NULL , showCurve=FALSE , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  #-----------------------------------------------------------------------------
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  theta = mcmcMat[,"theta"]
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  openGraph(width=4.0,height=3.0)
  par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  #-----------------------------------------------------------------------------
  postInfo = plotPost( theta , cex.lab = 1.75 , 
                       showCurve=showCurve ,
                       compVal=compVal , ROPE=rope , cex.main=1.5 ,
                       xlab=bquote(theta) , main=paste("theta") , 
                       col="skyblue" )
  z = sum(data$y)
  N = length(data$y)
  points( z/N , 0 , pch="+" , col="red" , cex=3 )
  text( max(theta) , 0 , bquote( z==.(z) ) , adj=c(1,-11) ) 
  text( max(theta) , 0 , bquote( N==.(N) ) , adj=c(1,-9.5) ) 
  
  #-----------------------------------------------------------------------------  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}

#===============================================================================
