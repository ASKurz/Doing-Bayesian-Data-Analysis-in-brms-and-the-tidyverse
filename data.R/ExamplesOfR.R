# ExamplesOfR.R
# Accompanies the book,
# Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
# A Tutorial with R, JAGS, and Stan. 2nd Ed. Academic Press.
# http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/

####################################################################
# Getting help:
?"?"  # help about help()
?"??" # help about help.search()

####################################################################
# Arithmetic operator precedence:
?"Syntax"     # information about operator precedence
1+2*3^2       # power first, then multiplication, then addition
(1+2)*3^2     # parentheses force addition before multiplication
(1+2*3)^2     # do operations inside parentheses before power
((1+2)*3)^2   # nested parentheses

####################################################################
# Logical operations:
?"!"
!TRUE  # negation
TRUE & FALSE  # conjunction
TRUE | FALSE  # disjunction
TRUE | TRUE & FALSE  # conjunction has precedence over disjunction
( TRUE | TRUE ) & FALSE  # parentheses force disjunction first

####################################################################
# Assignment, tests of equality
?"="
x <- 2  # assigns the value 2 to variable named x
x = 2   # assigns the value 2 to variable named x
x       # shows the contents of x
?"=="
x == 2  # checks whether the value of x is equal to 2
x != 3  # checks whether the value of x is NOT equal to 3
x < 3   # checks whether the value of x is less than 3
x > 3   # checks whether the value of x is greater than 3
x = 0.5 - 0.3
y = 0.3 - 0.1
x == y   # although mathematically TRUE, it's FALSE for limited precision
all.equal(x,y)  # equal up to precision of computer

####################################################################
# element-by-element vector operations:
?"c"
c(1,2,3) * c(7,6,5)
2 * c(1,2,3)
2 + c(1,2,3)

####################################################################
# colon operator and sequence function:
?":"
2+3:6      # colon operator has precedence over addition
(2+3):6    # parentheses override default precedence
1:3^2      # power operator has precedence over colon operator
(1:3)^2    # parentheses override default precedence
?"seq"
seq( from=0 , to=3 , by=0.5 )          # length not specified
seq( from=0 , to=3 , by=0.5001 )       # will not exceed end value
seq( from=0 , by=0.5 , length.out=7 )  # end not specified
seq( from=0 , to=3 , length.out=7 )    # increment not specified
seq( to=3 , by=0.5 , length.out=7 )    # start not specified

####################################################################
# replicate function:
?"rep"
ABC = c("A","B","C")  # define a vector for replication
rep( ABC, 2 )
rep( ABC, times=2 )
rep( ABC, times=c(4,2,1) )
rep( ABC, each=2 )
rep( ABC, each=2, length=10)
rep( ABC, each=2, times=3)
rep( ABC, each=2, times=c(1,2,3,1,2,3) )
rep( rep( ABC, each=2 ) , times=c(1,2,3,1,2,3) )

####################################################################
# getting at components of a vector:
x = c( 2.718 , 3.14 , 1.414 , 47405 )             # define the vector
names(x) = c( "e" , "pi" , "sqrt2" , "zipcode" )  # name the components
?"["
x[c(2,4)]                    # which indices to include
x[c(-1,-3)]                  # which indices to exclude
x[c(FALSE,TRUE,FALSE,TRUE)]  # for each position, include it?
x[c("pi","zipcode")]         # names of indices to include

####################################################################
# Factors:
?"factor"
x = c( "high" , "medium" , "low" , "high" , "medium" )
# Make it a factor, with default (alphabetical) levels:
xf = factor( x )
xf
as.numeric(xf)
# Re-order the levels:
xfo = factor( xf , levels=c("low","medium","high") , ordered=TRUE )
xfo
as.numeric(xfo)
# Labels:
xfol = factor( x , levels=c("low","medium","high") , ordered=TRUE ,
               labels=c("Bottom SES","Middle SES","Top SES") )
xfol

####################################################################
# Matrices:
?"matrix"
matrix( 1:6 , ncol=3 )  # contents are 1:6, filled by column
matrix( 1:6 , nrow=2 )  # or you can specify number of rows
matrix( 1:6 , nrow=2 , byrow=TRUE )  # filled by row instead of by column
matrix( 1:6 , nrow=2 ,       # with names of dimensions and rows and columns
        dimnames=list( TheRowDimName=c("Row1Name","Row2Name") ,
                       TheColDimName=c("Col1Name","Col2Name","Col3Name") ) )
x = matrix( 1:6 , nrow=2 ,   # assign matrix to x
            dimnames=list( TheRowDimName=c("Row1Name","Row2Name") ,
                           TheColDimName=c("Col1Name","Col2Name","Col3Name") ) )
x[2,3]  # use numerical indices
x["Row2Name","Col3Name"]  # use row, column names
x[2,1:3]  # specify range of columns for inclusion
x[2,]  # leave range of columns blank to include all columns
x[,3]  # all rows from column 3, returned as a vector
# The importance of the comma:
x[2,] # 2nd row (returned as vector)
x[,2] # 2nd column (returned as vector)
x[2] # no comma; returns 2nd element.

####################################################################
# Arrays:
?"array"
a = array( 1:24 , dim=c(3,4,2) , # 3 rows, 4 columns, 2 layers
           dimnames = list( RowDimName = c("R1","R2","R3") ,
                            ColDimName = c("C1","C2","C3","C4") ,
                            LayDimName = c("L1","L2") ) )
a
a["R3",,"L2"]  # returns all columns of R3 and L2, as a vector
a["R3","C4",]  # returns all layers of R3 and C4, as a vector

####################################################################
# Lists:
?"list"
MyList = list( "a"=1:3 , "b"=matrix(1:6,nrow=2) , "c"="Hello, world." )
MyList
MyList$a        # the contents of the list item named "a"
MyList$a[2]     # the second element of the list item named "a"
MyList[[1]]     # the contents of the first list item
MyList[[1]][2]  # the second element of the first list item
MyList[1]       # the first list item, including its name
MyList[1][2]    # does not make sense in this case

####################################################################
# Data frame:
?"data.frame"
d = data.frame( Integers=1:3 , NumberNames=c("one","two","three") )
d
d$NumberNames   # notice this is a factor
d[[2]]          # the second element contents
d[2]            # the second element with its name
d[,2]           # elements can be accessed as if it's a matrix
d[2,]           # elements can be accessed as if it's a matrix

####################################################################
# Reading data:
# HGN.csv contains the following:
#   Hair,Gender,Number,Name,Group
#   black,M,2,Alex,1
#   brown,F,4,Betty,1
#   blond,F,3,Carla,1
#   black,F,7,Diane,2
#   black,M,1,Edward,2
#   red,M,7,Frank,2
#   brown,F,10,Gabrielle,2
# The working directory must be set to 
# the folder in which that file is saved.
?"read.csv"
HGNdf = read.csv( "HGN.csv" )
HGNdf$Hair
as.numeric(HGNdf$Hair)
# Re-order the levels:
HGNdf$Hair = factor( HGNdf$Hair , levels=c("red","blond","brown","black") )
HGNdf$Hair
as.numeric(HGNdf$Hair)
# Change factor to vector:
HGNdf$Name
HGNdf$Name = as.vector( HGNdf$Name )
HGNdf$Name
# Change vector to factor:
HGNdf$Group
HGNdf$Group = factor( HGNdf$Group )
HGNdf$Group

####################################################################
# Saving data:
?"write.csv"
write.csv( HGNdf , file="HGN.csv" , row.names=FALSE , quote=FALSE )
?"save"
save( HGNdf , file="HGN.Rdata" )
?"load"
load( "HGN.Rdata" )
?"objects"
objects()

####################################################################
# Utility functions:

?"summary"
x = c( rep(1,100) , rep(2,200) , rep(3,300) )  # 100 1's, 200 2's, 300 3's
summary(x)
xf = factor(x)
summary(xf)

?"head"
?"tail"
?"str"

?"aggregate"
aggregate( x=HGNdf$Number , by=list(HGNdf$Gender,HGNdf$Hair) , FUN=median )
# explicitly name the variables:
aggregate( x=list(Number=HGNdf$Number) ,
           by=list(Gender=HGNdf$Gender,Hair=HGNdf$Hair) , FUN=median )
# use formula format:
aggregate( Number ~ Gender + Hair , data=HGNdf , FUN=median )

# Table of counts:
# Using aggregate, to produce one count per row with factor levels each row:
aggregate( x=list(Count=rep(1,NROW(HGNdf))) , # column of 1's
           by=list(Gender=HGNdf$Gender,Hair=HGNdf$Hair) , FUN=sum )
# Using table, to produce several counts per row with implicit factor levels:
?"table"
table(list(Gender=HGNdf$Gender,Hair=HGNdf$Hair))

?"apply"
# Display the array, a, defined above:
a
# Sum within columns and layers, collapsed across rows, i.e.,
# retain the 2nd and 3rd dimensions, while collapsing across the other dimension. 
apply( a , MARGIN=c(2,3) , FUN=sum )

####################################################################
# Re-shaping a data arrangement
install.packages("reshape2") # only needed once
library(reshape2)

# Display the array, a, defined above:
a
# Re-shape array as one datum per row, with dimension-level identifiers:
?melt
am = melt(a) # result is data frame
am

####################################################################
# Running a saved program.
# The file SimpleGraph.R contains these lines:
#   x = seq( from = -2 , to = 2 , by = 0.1 )
#   y = x^2                                 
#   plot( x , y , type = "l" )              
# The working directory must be set to 
# the folder in which that file is saved.
?"source"
source("SimpleGraph.R")

####################################################################
# User-defined functions:
?"function"
asqplusb = function( a , b=1 ) {  # begin function definition
  c = a^2 + b
  return( c )
}  # end of function definition
# try the function:
asqplusb( a=3 , b=2 )
# explicit argument labels may go in any order:
asqplusb( b=2 , a=3 )
# implicit argument labels must be in defined order:
asqplusb( 3 , 2 )
asqplusb( 2 , 3 )
# arguments with defined default values need not be specified:
asqplusb( a=2 , b=1 )
asqplusb( a=2 )  # b gets default value
asqplusb( b=1 )  # error: a has no default value
asqplusb( 2 )    # value is assigned to first argument

####################################################################
# Conditions and loops
?"Control"
x = 5
if ( x <= 3 ) {  # if x is less than or equal to 3
  show("small")  # display the word "small"
} else {         # otherwise
  show("big")    # display the word "big"
}                # end of 'else' clause
# The placement of curly brace before 'else' matters:
if ( x <= 3 ) { show("small") }
else { show("big") }  # error!

for ( countDown in 5:1 ) {
  show(countDown)
}
for ( note in c("do","re","mi") ) {
  show(note)
}

?"Control"

####################################################################
# Measuring processing time

startTime = proc.time()
y = vector(mode="numeric",length=1.0E6)
for ( i in 1:1.0E6 ) { y[i] = log(i) }
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

startTime = proc.time()
y = log(1:1.0E6)
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

####################################################################
# opening graphics windows and saving graphs
source("DBDA2E-utilities.R")             # read defn. of openGraph, saveGraph
openGraph( width=3 , height=4 )          # open a graphics window
plot( x=1:4 , y=c(1,3,2,4) , type="o" )  # make a plot in the screen window
saveGraph( file="temp" , type="pdf" )    # save the graph as "temp.pdf"
