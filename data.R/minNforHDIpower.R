minNforHDIpower = function( genPriorMode , genPriorN ,
                            HDImaxwid=NULL , nullVal=NULL , 
                            ROPE=c(max(0,nullVal-.02),min(1,nullVal+.02)) ,
                            desiredPower=0.8 , audPriorMode=0.5 , audPriorN=2 ,
                            HDImass=0.95 , initSampSize=20 , verbose=TRUE ) {
  # Check for argument consistency:
  if ( !xor( is.null(HDImaxwid) , is.null(nullVal) ) ) {
    stop("One and only one of HDImaxwid and nullVal must be specified.")
  }
  # Load HDIofICDF function if not already present:
  if ( !exists("HDIofICDF") ) source("DBDA2E-utilities.R")
  # Convert prior mode and N to a,b parameters of beta distribution:
  genPriorA = genPriorMode * (genPriorN-2) + 1
  genPriorB = ( 1.0 - genPriorMode ) * (genPriorN-2) + 1
  audPriorA = audPriorMode * (audPriorN-2) + 1
  audPriorB = ( 1.0 - audPriorMode ) * (audPriorN-2) + 1
  # Initialize loop for incrementing sampleSize:
  sampleSize = initSampSize
  notPowerfulEnough = TRUE
  # Increment sampleSize until desired power is achieved:
  while( notPowerfulEnough ) {
    zvec = 0:sampleSize # vector of all possible z values for N flips.
    # Compute probability of each z value for data-generating prior:
    pzvec = exp( lchoose( sampleSize , zvec )
                 + lbeta( zvec + genPriorA , sampleSize-zvec + genPriorB )
                 - lbeta( genPriorA , genPriorB ) )
    # For each z value, compute posterior HDI: 
    # hdiMat will hold HDI limits for each z:
    hdiMat = matrix( 0 , nrow=length(zvec) , ncol=2 )
    for ( zIdx in 1:length(zvec) ) {
      z = zvec[zIdx]
      hdiMat[zIdx,] = HDIofICDF( qbeta , 
                                 shape1 = z + audPriorA ,
                                 shape2 = sampleSize - z + audPriorB ,
                                 credMass = HDImass )
    }
    # Compute HDI widths:
    hdiWid = hdiMat[,2] - hdiMat[,1]
    # Sum the probabilities of outcomes with satisfactory HDI widths:
    if ( !is.null( HDImaxwid ) ) {
      powerHDI = sum( pzvec[ hdiWid < HDImaxwid ] )
    }
    # Sum the probabilities of outcomes with HDI excluding ROPE:
    if ( !is.null( nullVal ) ) {
      powerHDI = sum( pzvec[ hdiMat[,1] > ROPE[2] | hdiMat[,2] < ROPE[1] ] )
    }
    if ( verbose ) {
      cat( " For sample size = ", sampleSize , ", power = " , powerHDI ,
           "\n" , sep="" ) ; flush.console()
    }
    if ( powerHDI > desiredPower ) {  # If desired power is attained,
      notPowerfulEnough = FALSE       # set flag to stop,
    } else {                          # otherwise
      sampleSize = sampleSize + 1     # increment the sample size.
    }
  } # End while( notPowerfulEnough ).
  # Return the sample size that achieved the desired power:
  return( sampleSize )
} # End of function.
