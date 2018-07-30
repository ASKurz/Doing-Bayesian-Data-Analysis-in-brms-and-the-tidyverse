# Jags-Ydich-Xnom1subj-MbernBetaModelComp.R
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
graphics.off()
rm(list=ls(all=TRUE))
source("DBDA2E-utilities.R")
require(rjags)
fileNameRoot="Jags-Ydich-Xnom1subj-MbernBetaModelCompPseudoPrior-PSEUDO-" 

#------------------------------------------------------------------------------
# THE DATA.

N=30
z=ceiling(.55*N)
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
  theta <- equals(m,1)*theta1 + equals(m,2)*theta2
  theta1 ~ dbeta( omega1[m]*(kappa1[m]-2)+1 , (1-omega1[m])*(kappa1[m]-2)+1 ) 
  omega1[1] <- .10 # true prior value
  omega1[2] <- .40 # pseudo prior value
  kappa1[1] <- 20 # true prior value
  kappa1[2] <- 50 # pseudo prior value
  theta2 ~ dbeta( omega2[m]*(kappa2[m]-2)+1 , (1-omega2[m])*(kappa2[m]-2)+1 ) 
  omega2[1] <- .70 # pseudo prior value
  omega2[2] <- .90 # true prior value
  kappa2[1] <- 50 # pseudo prior value
  kappa2[2] <- 20 # true prior value
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
# initsList = list( theta1=0.5 , theta2=0.5 , m=1 ) 

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("m","theta1","theta2") 
adaptSteps = 1000            # Number of steps to "tune" the samplers.
burnInSteps = 1000           # Number of steps to "burn-in" the samplers.
nChains = 4                  # Number of chains to run.
numSavedSteps=10000          # Total number of steps in chains to save.
thinSteps=1                  # Number of steps to "thin" (1=keep every step).
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
for ( parName in c("m") ) { # parameterNames ) {
  diagMCMC( codaSamples , parName=parName ,
            saveName=fileNameRoot , saveType="eps" )
}

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcMat = as.matrix( codaSamples , chains=TRUE )
m = mcmcMat[,"m"]
# Compute the proportion of m at each index value:
pM1 = sum( m == 1 ) / length( m )
pM2 = 1 - pM1
# Extract theta values for each model index:
theta1M1 = mcmcMat[,"theta1"][ m == 1 ] # true theta1
theta1M2 = mcmcMat[,"theta1"][ m == 2 ] # pseudo theta1
theta2M1 = mcmcMat[,"theta2"][ m == 1 ] # pseudo theta2
theta2M2 = mcmcMat[,"theta2"][ m == 2 ] # true theta2

# Plot histograms of sampled theta values for each model,
# with pM displayed.
openGraph(width=7,height=7)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,1,2,3,4,5),nrow=3,byrow=TRUE)   )
plotPost( m , breaks=seq(0.95,2.05,0.1) , xlim=c(0.75,2.25) , 
          cenTend="mean" , xlab="m" , cex.main=1.75 ,
          main=bquote( "Model Index." * 
                         " p(m=1|D) =" * .(signif(pM1,3)) *
                         ", p(m=2|D) =" * .(signif(pM2,3)) ) )
plotPost( theta1M1 , 
          main=bquote( theta[1]*" when m=1 (using true prior)" ) , 
          cex.main=1.75 , xlab=bquote(theta[1]) , xlim=c(0,1) )
plotPost( theta2M1 , 
          main=bquote( theta[2]*" when m=1; pseudo-prior" ) , 
          cex.main=1.75 , xlab=bquote(theta[2]) , xlim=c(0,1) )
plotPost( theta1M2 , 
          main=bquote( theta[1]*" when m=2; pseudo-prior" ) , 
          cex.main=1.75 , xlab=bquote(theta[1]) , xlim=c(0,1) )
plotPost( theta2M2 , 
          main=bquote( theta[2]*" when m=2 (using true prior)" ) , 
          cex.main=1.75 , xlab=bquote(theta[2]) , xlim=c(0,1) )
saveGraph( file=paste0(fileNameRoot,"Post") , type="eps" )
