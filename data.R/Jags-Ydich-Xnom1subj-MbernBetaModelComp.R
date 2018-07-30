# Jags-Ydich-Xnom1subj-MbernBetaModelComp.R
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
graphics.off()
rm(list=ls(all=TRUE))
source("DBDA2E-utilities.R")
require(rjags)
fileNameRoot="Jags-Ydich-Xnom1subj-MbernBetaModelComp-" # for output filenames

#------------------------------------------------------------------------------
# THE DATA.

N=9
z=6
y = c( rep(0,N-z) , rep(1,z) )
dataList = list(
  y = y ,
  N = N 
)

#------------------------------------------------------------------------------
# THE MODEL.

modelString = "
model {
  for ( i in 1:N ) {
    y[i] ~ dbern( theta )
  }
  theta ~ dbeta( omega[m]*(kappa-2)+1 , (1-omega[m])*(kappa-2)+1 ) 
  omega[1] <- .25
  omega[2] <- .75
  kappa <- 12
  m ~ dcat( mPriorProb[] )
  mPriorProb[1] <- .5
  mPriorProb[2] <- .5
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Specific initialization is not necessary in this case, 
# but here is a lazy version if wanted:
# initsList = list( theta=0.5 , m=1 ) 

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("theta","m") 
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 1000           # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
numSavedSteps=50000          # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , # inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

#------------------------------------------------------------------------------- 
# Display diagnostics of chain:

parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaSamples , parName=parName ,
            saveName=fileNameRoot , saveType="eps" )
}

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcMat = as.matrix( codaSamples , chains=TRUE )
m = mcmcMat[,"m"]
theta = mcmcMat[,"theta"]

# Compute the proportion of m at each index value:
pM1 = sum( m == 1 ) / length( m )
pM2 = 1 - pM1

# Extract theta values for each model index:
thetaM1 = theta[ m == 1 ]
thetaM2 = theta[ m == 2 ]

# Plot histograms of sampled theta values for each model,
# with pM displayed.
openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,1,2,3),nrow=2,byrow=FALSE) , widths=c(1,2) )
plotPost( m , breaks=seq(0.9,2.1,0.2) , cenTend="mean" , xlab="m" , main="Model Index" )
plotPost( thetaM1 , 
          main=bquote( theta*" when m=1" * " ; p(m=1|D)" == .(signif(pM1,3)) ) , 
          cex.main=1.75 , xlab=bquote(theta) , xlim=c(0,1) )
plotPost( thetaM2 , 
          main=bquote( theta*" when m=2" * " ; p(m=2|D)" == .(signif(pM2,3)) ) , 
          cex.main=1.75 , xlab=bquote(theta) , xlim=c(0,1) )
saveGraph( file=paste0(fileNameRoot,"Post") , type="eps" )
