 
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

fileNameRoot="Stan-BernBeta-Script-"

source("DBDA2E-utilities.R")
library(rstan)

# Specify model:
modelString = "
  data {
    int<lower=0> N ;
    int y[N] ; 
  }
  parameters {
    real<lower=0,upper=1> theta ;
  }
  model {
    theta ~ beta(1,1) ;
    y ~ bernoulli(theta) ; 
  }
" # close quote for modelString

# Translate model to C++ and compile to DSO:
stanDso <- stan_model( model_code=modelString ) 

# Specify data:
N = 50 ; z = 10
y = c(rep(1,z),rep(0,N-z))
dataList = list(
  y = y ,
  N = N 
)

# Generate posterior sample:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     chains = 3 ,
                     iter = 1000 , 
                     warmup = 200 , 
                     thin = 1 )

openGraph()
traceplot(stanFit,pars=c("theta"))
saveGraph(file=paste0(fileNameRoot,"StanTrace"),type="eps")
openGraph()
plot(stanFit,pars=c("theta"))
saveGraph(file=paste0(fileNameRoot,"StanPlot"),type="eps")

# Make graphs:
# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format. This excludes warmup and thinned steps.
mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , 
                              function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
diagMCMC( mcmcCoda , parName=c("theta") )
saveGraph(file=paste0(fileNameRoot,"Diag"),type="eps")

#------------------------------------------------------------------------------
# Another data set:

# Specify data:
N = 50 ; z = 40
y = c(rep(1,z),rep(0,N-z))
dataList = list(
  y = y ,
  N = N 
)

# Generate posterior sample:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     chains = 3 ,
                     iter = 1000 , 
                     warmup = 200 , 
                     thin = 1 )

openGraph()
traceplot(stanFit,pars=c("theta"))
openGraph()
plot(stanFit,pars=c("theta"))

# Make graphs:
# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format. This excludes warmup and thinned steps.
mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , 
                              function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
diagMCMC( mcmcCoda , parName=c("theta") )

#------------------------------------------------------------------------------
# Generate sample from prior:

# Specify model:
modelString = "
  data {
    int<lower=0> N ;
    int y[N] ; 
  }
  parameters {
    real<lower=0,upper=1> theta ;
  }
  model {
    theta ~ beta(1,1) ;
//    y ~ bernoulli(theta) ;  // likelihood commmented out
  }
" # close quote for modelString

# Translate model to C++ and compile to DSO:
stanDso <- stan_model( model_code=modelString ) 

# Specify data:
N = 50 ; z = 10
y = c(rep(1,z),rep(0,N-z))
dataList = list(
  y = y ,
  N = N 
)

# Generate posterior sample:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     chains = 3 ,
                     iter = 1000 , 
                     warmup = 200 , 
                     thin = 1 )

openGraph()
traceplot(stanFit,pars=c("theta"))
openGraph()
plot(stanFit,pars=c("theta"))

# Make graphs:
# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format. This excludes warmup and thinned steps.
mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , 
                              function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
diagMCMC( mcmcCoda , parName=c("theta") )
