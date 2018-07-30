# Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa-Power.R
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
set.seed(47405) # Optional, merely for replicability.

# Load the functions genMCMC, smryMCMC, and plotMCMC:
# (This also sources DBDA2E-utilities.R)
source("Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa.R")

# Specify idealized hypothesis:
idealGroupMean = 0.65
idealGroupSD = 0.07
idealNsubj = 100       # more subjects => higher confidence in hypothesis
idealNtrlPerSubj = 100 # more trials => higher confidence in hypothesis
# Generate random theta values for idealized subjects:
betaAB = betaABfromMeanSD( idealGroupMean , idealGroupSD )
theta = rbeta( idealNsubj , betaAB$a , betaAB$b )
# Transform the theta values to exactly match idealized mean, SD:
theta = ((theta-mean(theta))/sd(theta))*idealGroupSD + idealGroupMean
theta[ theta >= 0.999 ] = 0.999 # must be between 0 and 1
theta[ theta <= 0.001 ] = 0.001 # must be between 0 and 1
# Generate idealized data very close to theta's:
z = round( theta*idealNtrlPerSubj ) 
# Convert to data format needed by JAGS function:
dataMat=matrix(0,ncol=2,nrow=0,dimnames=list(NULL,c("y","s")))
for ( sIdx in 1:idealNsubj ) {
  yVec = c(rep(1,z[sIdx]),rep(0,idealNtrlPerSubj-z[sIdx]))
  dataMat = rbind( dataMat , cbind( yVec , rep(sIdx,idealNtrlPerSubj) ) )
}
idealDatFrm = data.frame(dataMat)
# Run Bayesian analysis on idealized data:
mcmcCoda = genMCMC( data=idealDatFrm , saveName=NULL , 
                    numSavedSteps=2000 , thinSteps=20 )
# Convert coda object to matrix for convenience:
mcmcMat = as.matrix(mcmcCoda)

# # Examine idealized parameter values:
# effectiveSize(mcmcCoda[,"omega"])
# effectiveSize(mcmcCoda[,"kappa"])
openGraph(width=7,height=3)
layout(matrix(1:2,ncol=2))
par( mar=c(3,1,1,1) , mgp=c(2.0,0.7,0) ,  oma=0.1+c(0,0,2,0) ,
     cex.lab=1.75 , cex.main=1.5 ,  pch=20 )
plotPost( mcmcMat[,"omega"] , xlab="omega" , cenTend="mean" , xlim=c(.45,.85) ,
          border="skyblue" , HDItextPlace=0.9 )
plotPost( mcmcMat[,"kappa"] , xlab="kappa" , cenTend="mode" , xlim=c(0,250) )
mtext( text=bquote(list( idealNsubj==.(idealNsubj) , 
                         idealNtrlPerSubj==.(idealNtrlPerSubj)  )) , 
       outer=TRUE , adj=c(0.5,0.5) , cex=1.5 )
saveGraph( file=paste0("Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa-Power-",
                       idealNsubj,"-",idealNtrlPerSubj) , type="eps" )

# Define function that assays goal achievement for a single set of data:
goalAchievedForSample = function( data ) {
  # Generate the MCMC chain:
  mcmcCoda = genMCMC( data=data , saveName=NULL ,
                      numSavedSteps=5000 , thinSteps=2 ) 
                      #numSavedSteps=11000 , thinSteps=20  
  # Convert coda object to matrix for convenience:
  mcmcMat = as.matrix(mcmcCoda)
  # Specify criteria for goals:
  nullROPE = c(0.48,0.52)
  HDImaxWid = 0.2  
  # Compute HDIs:
  HDImat = apply( mcmcMat , 2 , "HDIofMCMC" )
  show( HDImat[,1:5] )
  # Define list for recording results:
  goalAchieved = list()
  # Goal: omega greater than ROPE:
  goalAchieved = c( goalAchieved , 
                    "omegaAboveROPE"=unname( HDImat[1,"omega"] > nullROPE[2] ) )
  # Goal: omega HDI width less than max width:
  goalAchieved = c( goalAchieved , 
                    "omegaNarrowHDI"=unname( HDImat[2,"omega"]-HDImat[1,"omega"] 
                                             < HDImaxWid ) )
  # Goal: at least one theta greater than ROPE with none below:
  thetaCols = grep("theta",colnames(HDImat)) # column indices of thetas
  goalAchieved = c( goalAchieved , 
                    "thetasAboveROPE"= (any(HDImat[1,thetaCols] > nullROPE[2]) 
                                        & !any(HDImat[2,thetaCols] < nullROPE[1])))
  # Goal: all theta's HDI width less than max width:
  goalAchieved = c( goalAchieved , 
                    "thetasNarrowHDI"= all( HDImat[2,thetaCols] 
                                            - HDImat[1,thetaCols] 
                                            < HDImaxWid ) )
  # More goals can be inserted here if wanted...
  # Return list of goal results:
  return(goalAchieved)
}

# Specify sample size for each simulated data set:
Nsubj = 2*7 ; NtrlPerSubj = 47  # 658 flips total
#Nsubj = 7 ; NtrlPerSubj = 2*47  # 658 flips total
# Specify the number of simulated experiments:
nSimulatedDataSets = min(500,NROW(mcmcMat)) # An arbitrary large number.
# Run the simulated experiments:
simCount=0
if (exists("goalTally")) rm(goalTally) # in case previously run from here down
for ( simIdx in ceiling(seq(1,NROW(mcmcMat),length=nSimulatedDataSets)) ) {
  simCount=simCount+1
  cat( "\n\n==================== Simulation",simCount,"of",nSimulatedDataSets,
       "====================\n\n" ) 
  # Generate random omega and kappa for group distribution:
  genOmega = mcmcMat[simIdx,"omega"]
  genKappa = mcmcMat[simIdx,"kappa"]
  # Generate random theta's for individuals:
  genTheta = rbeta( Nsubj , genOmega*(genKappa-2)+1 , (1-genOmega)*(genKappa-2)+1 )
  # Generate random data based on parameter value:
  dataMat=matrix(0,ncol=2,nrow=0,dimnames=list(NULL,c("y","s")))
  for ( sIdx in 1:Nsubj ) {
    z = rbinom( 1 , size=NtrlPerSubj , prob=genTheta[sIdx] )
    yVec = c(rep(1,z),rep(0,NtrlPerSubj-z))
    dataMat = rbind( dataMat , cbind( yVec , rep(sIdx,NtrlPerSubj) ) )
  }
  # Do Bayesian analysis on simulated data:
  goalAchieved = goalAchievedForSample( data.frame(dataMat) )
  # Tally the results:
  if (!exists("goalTally")) { # if goalTally does not exist, create it
    goalTally=matrix( nrow=0 , ncol=length(goalAchieved) ) 
  }
  goalTally = rbind( goalTally , goalAchieved )
  # save( goalTally ,
  #       file="Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa-Power-goalTally.Rdata" )
}

# For each goal...
for ( goalIdx in 1:NCOL(goalTally) ) {
  # Extract the goal name for subsequent display:
  goalName = colnames(goalTally)[goalIdx]
  # Compute number of successes:
  goalHits = sum(unlist(goalTally[,goalIdx]))
  # Compute number of attempts:
  goalAttempts = NROW(goalTally)
  # Compute proportion of successes:
  goalEst = goalHits/goalAttempts
  # Compute HDI around proportion:
  goalEstHDI = HDIofICDF( qbeta ,
                          shape1=1+goalHits , 
                          shape2=1+goalAttempts-goalHits )
  # Display the result:
  show( paste0( goalName,
                ": Est.Power=" , round(goalEst,3) , 
                "; Low Bound=" , round(goalEstHDI[1],3) ,
                "; High Bound=" , round(goalEstHDI[2],3) ) )
}
