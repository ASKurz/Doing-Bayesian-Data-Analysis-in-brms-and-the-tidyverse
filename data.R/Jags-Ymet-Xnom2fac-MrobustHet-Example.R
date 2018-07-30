# Example for Jags-Ymet-Xnom2fac-MrobustHet.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#Load The data file 

fileNameRoot = "SalaryRobustHet-" 
graphFileType = "eps" 
myDataFrame = read.csv( file="Salary.csv" )
# Re-label and re-order the Pos factor:
myDataFrame$Pos = factor( myDataFrame$Pos , 
                          levels=c("FT3","FT2","FT1","NDW","DST") , 
                          ordered=TRUE , 
                          labels=c("Assis","Assoc","Full","Endow","Disting") )
# Specify the column names in the data file relevant to the analysis:
yName="Salary" 
# x1 should be factor with fewer levels, to plot in single pane:
x1Name="Pos" 
x2Name="Org" 
# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):

x1contrasts = list( 
  list( c("Full") , c("Assoc") , compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( c("Assoc") , c("Assis") , compVal=0.0 , ROPE=c(-1000,1000) ) 
)
x2contrasts = list( 
  list( c("CHEM") , c("ENG") , compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( c("CHEM") , c("PSY") , compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( c("BFIN") , c("PSY","CHEM","ENG") , compVal=0.0 , ROPE=c(-1000,1000) ) 
)
# Each interaction contrast is a list of 2 lists of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL)::
x1x2contrasts = list( 
  list( list( c("Full") , c("Assis") ) ,
        list( c("CHEM") , c("ENG") ) ,
        compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( list( c("Full") , c("Assis") ) ,
        list( c("CHEM") , c("PSY") ) ,
        compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( list( c("Full") , c("Assoc","Assis") ) ,
        list( c("BFIN") , c("PSY","CHEM","ENG") ) , 
        compVal=0.0 , ROPE=c(-1000,1000) )
) 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2fac-MrobustHet.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=15000 , thinSteps=5 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("b0","b1[1]","b2[1]","b1b2[1,1]","ySigma[1,1]","ySigma[1,7]","ySigma[5,7]","nu") ) {
    diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , x1Name=x1Name , x2Name=x2Name ,
                        x1contrasts=x1contrasts , 
                        x2contrasts=x2contrasts , 
                        x1x2contrasts=x1x2contrasts ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name ,
          x1contrasts=x1contrasts , 
          x2contrasts=x2contrasts , 
          x1x2contrasts=x1x2contrasts ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
# Other specific comparisons of cells:
if ( fileNameRoot == "SalaryRobustHet-" ) {
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Full"
  THATx1 = "Assis"
  ATx2 = "CHEM"
  THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
  THATidx = which(levels(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Full"
  THATx1 = "Assis"
  ATx2 = "PSY"
  THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
  THATidx = which(levels(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x2level minus THAT x2level at AT x1level:
  THISx2 = "PSY"
  THATx2 = "ENG"
  ATx1 = "Full"
  THISidx = which(levels(myDataFrame[,x2Name])==THISx2)
  THATidx = which(levels(myDataFrame[,x2Name])==THATx2)
  ATidx   = which(levels(myDataFrame[,x1Name])==ATx1)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THISidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THATidx,"]",sep="")] , 
    main=paste(THISx2,"-",THATx2,"@",ATx1) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx2,"-",THATx2,"At",ATx1,sep=""),
            type=graphFileType)
}
#------------------------------------------------------------------------------- 
