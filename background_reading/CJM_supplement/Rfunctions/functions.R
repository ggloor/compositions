# rdirichlet function used in ALDEx2
# extracted here since it is not exposed
#copied from mc2d R package
#licenced GPL>=2
#should be compatable with AGPL3
rdirichlet <- function (n, alpha)
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  #if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  #n <- as.integer(n)
  if(n < 0) stop("value(n) can not be negative in rtriang")

  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}

#######################################################################################
#' @title Expected value of phi from Dirichlet log-ratio distributions
#' Copyright Greg Gloor, 2016
#' Licensed under AGPL3 license
#' @description
#' returns dataframe of the lower-triangle of symmetrical phi metric
#' where the value of phi is the expected value of a number of Dirichlet
#' Monte-Carlo replicates of the data. This reduces the problem of
#' 0-count and low-count features being highly variable because their
#' values range wildly and so the expected value is always large
#' @details requires aldex.clr function from ALDEx2
#' param aldex.clr is an S3 object from the aldex.clr function
#' we ignore all the other measures that are used for trouble-shooting phi
#' the sma.df function in particular is very time and memory intensive
#' @examples
#' # use a count table where the samples are by column, features by row
#' x <- aldex.clr(count.table)
#' # returns a dataframe of the expected value of the lower triangle of the
#' propr.phisym function. The number of Dirichlet Monte-Carlo replicates is
#' obtained from the aldex.clr object

aldex.phi <- function(aldex.clr){

	# calculate expected value of phi
	# a single high phi value will push the component out of consideration
	# a median is right out for memory considerations

	# get first value
	sym.phi <- propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
	    function(y){y[,1]})))

	# sum the rest of the values as we proceed through the DIR MC instances
	for(i in 2:numMCInstances(aldex.clr)){
		#print(i)
		sym.phi <- sym.phi + propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
		    function(y){y[,i]})))
	}
	##### Done ALDEx2 stuff

	# make indices of the correct size
	lt <- which(col(sym.phi)<row(sym.phi), arr.ind=FALSE)
	lt.ind <- which(col(sym.phi)<row(sym.phi), arr.ind=TRUE)

	# dataframe to hold the info,
	# data is a set of vectors where only the lower triangle is kept, because the matrix
	#    is symmetrical
	# this is needed so subsequent subset function works properly
	sma.df <- data.frame(row=factor(rownames(sym.phi)[lt.ind[,"row"]]),
		col=factor(colnames(sym.phi)[lt.ind[,"col"]]))

	#save the lower triangle as an expected value
	sma.df$phi <- sym.phi[lt] /  numMCInstances(aldex.clr)

	return(sma.df)
}
