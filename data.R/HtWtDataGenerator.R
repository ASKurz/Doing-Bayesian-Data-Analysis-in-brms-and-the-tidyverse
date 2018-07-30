HtWtDataGenerator = function( nSubj , rndsd=NULL , maleProb=0.50 ) {
  # Random height, weight generator for males and females. Uses parameters from
  # Brainard, J. & Burmaster, D. E. (1992). Bivariate distributions for height and
  # weight of men and women in the United States. Risk Analysis, 12(2), 267-275.
  # Kruschke, J. K. (2011). Doing Bayesian data analysis:
  # A Tutorial with R and BUGS. Academic Press / Elsevier.
  # Kruschke, J. K. (2014). Doing Bayesian data analysis, 2nd Edition:
  # A Tutorial with R, JAGS and Stan. Academic Press / Elsevier.
  
  require(MASS)
  
  # Specify parameters of multivariate normal (MVN) distributions.
  # Men:
  HtMmu = 69.18
  HtMsd = 2.87
  lnWtMmu = 5.14
  lnWtMsd = 0.17
  Mrho = 0.42
  Mmean = c( HtMmu , lnWtMmu )
  Msigma = matrix( c( HtMsd^2 , Mrho * HtMsd * lnWtMsd ,
                      Mrho * HtMsd * lnWtMsd , lnWtMsd^2 ) , nrow=2 )
  # Women cluster 1:
  HtFmu1 = 63.11
  HtFsd1 = 2.76
  lnWtFmu1 = 5.06
  lnWtFsd1 = 0.24
  Frho1 = 0.41
  prop1 = 0.46
  Fmean1 = c( HtFmu1 , lnWtFmu1 )
  Fsigma1 = matrix( c( HtFsd1^2 , Frho1 * HtFsd1 * lnWtFsd1 ,
                       Frho1 * HtFsd1 * lnWtFsd1 , lnWtFsd1^2 ) , nrow=2 )
  # Women cluster 2:
  HtFmu2 = 64.36
  HtFsd2 = 2.49
  lnWtFmu2 = 4.86
  lnWtFsd2 = 0.14
  Frho2 = 0.44
  prop2 = 1 - prop1
  Fmean2 = c( HtFmu2 , lnWtFmu2 )
  Fsigma2 = matrix( c( HtFsd2^2 , Frho2 * HtFsd2 * lnWtFsd2 ,
                       Frho2 * HtFsd2 * lnWtFsd2 , lnWtFsd2^2 ) , nrow=2 )
  
  # Randomly generate data values from those MVN distributions.
  if ( !is.null( rndsd ) ) { set.seed( rndsd ) }
  datamatrix = matrix( 0 , nrow=nSubj , ncol=3 )
  colnames(datamatrix) = c( "male" , "height" , "weight" )
  maleval = 1 ; femaleval = 0 # arbitrary coding values
  for ( i in 1:nSubj ) {
    # Flip coin to decide sex
    sex = sample( c(maleval,femaleval) , size=1 , replace=TRUE , 
                  prob=c(maleProb,1-maleProb) )
    if ( sex == maleval ) { datum = mvrnorm(n = 1, mu=Mmean, Sigma=Msigma ) }
    if ( sex == femaleval ) {
      Fclust = sample( c(1,2) , size=1 , replace=TRUE , prob=c(prop1,prop2) )
      if ( Fclust == 1 ) { datum = mvrnorm(n = 1, mu=Fmean1, Sigma=Fsigma1 ) }
      if ( Fclust == 2 ) { datum = mvrnorm(n = 1, mu=Fmean2, Sigma=Fsigma2 ) }
    }
    datamatrix[ i , ] = c( sex , round( c( datum[1] , exp( datum[2] ) ) , 1 ) )
  }
  
  return( datamatrix )
} # end function