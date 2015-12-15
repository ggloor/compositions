# @title Filter compositional dataset
#
# @ description
# \code{codamb.filter} returns a reduced able of counts where the
# samples must contain a minimum number of reads, and OTUs must be
# found with a minimum abundance in all remaining samples.
# \emph{Note:} filters by min.prop first, min.fraction second, min.prop third.
# Requires numeric data only.
# \emph{Note:} There is a parallel function (to be implemented) that
# reduces the metadata with the same filtering parameters
#
# @param x A matrix or dataframe containing a count table
# @param min.reads The minimum reads per sample. Default=5000
# @param min.prop The minimum proportion of a read in any sample. Default 0.001
# @param min.fraction The minimum sample proportion of non-0 reads for each variable
# @param samples.by.row True if rows contain samples, False if rows contain variables
#
# @return A reduced set of data with samples by rows
#
# @examples
# to do
# @ export
codamb.filter <- function(x, y=tax.vector, min.reads=5000, min.prop=0.001,
  min.fraction=0, samples.by.row=TRUE){
  if(samples.by.row==FALSE) data <-x
  if(samples.by.row==TRUE) data <- t(x)
  # todo: check for numeric
  data.0 <- data[,which(apply(data,2,sum) > min.reads)]

  d.frac <- apply(data.0, 2, function(x){x/sum(x)})
  data.1 <- data.0[which(apply(d.frac, 1, max) > min.prop),]
  rm(d.frac)

  return( data.frame(data.1[which(apply(data.1, 1,
    function(x){length(which(x != 0))/length(x)}) > min.fraction),]) )
}

# @title Center log-ratio function
#
# @description
# \code{codamb.clr} returns a matrix of center log-ratio transformed data
# with samples by row
# equivalent to log(x/gx) for every value where gx is the geomtric mean of the vector X
# \emph{Note:} natural log is used for biplots and other exploratory analyses
# @param x A matrix or dataframe with samples by row or column
# @param samples.by.row TRUE if samples by row, FALSE if samples by column
# @return Center log-ratio transform of the data
codamb.clr <- function(x, samples.by.row=TRUE){
  if(min(x) == 0) stop("0 values must be replaced, estimated, or eliminated")
  if(samples.by.row == TRUE) margin=1
  if(samples.by.row == FALSE) margin=2

  return ( t(apply(d.n0, margin, function(x){log(x) - mean(log(x))})) )
}

# @title Identifying sample outliers
#
# @description
# \code{codamb.outlier} returns a vector of proportional contribution to group variance
# \emph{Note:} Samples must be grouped. This approach makes no sense across
# groups. If you do not know if you have natural groups, ignore this step and
# examine your data by PCA
# \emph{Note:} Requires the \code{compositions} R package
# @param x A matrix or dataframe with clr transformed values, with samples by row
# @return vector of proportional variance contributions for each sample
codamb.outlier <- function(x){
  mv <- mvar(x)
  pcx <- prcomp(x)

  return ( apply(pcx$x,1,function(y){sum(y^2/mv)}) )
}



