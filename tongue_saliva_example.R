
############### EXAMPLE CODE
# required libraries for compositional analysis
library(compositions)
library(zCompositions)

# read the entire data table
# tongue samples start with td_
# saliva samples start with sa_
d <- read.table("tongue_saliva.txt", header=T, row.names=1, sep="\t")

# get the taxon list and remove that column from the table
tax <- d$tax.0
d$tax.0 <- NULL

# samples must be by rows so use t()
d.n0 <- cmultRepl(t(d), label=0, method="CZM")

# the apply function rotates the data when by row, so turn it back to samples by row
# how I hate R
d.n0.clr <- t(apply(d.n0, 1, function(x){log2(x) - mean(log2(x))}))

####
# this set of code replaces row and column names with single values
# and the genus name of each taxon, if known
# replace rownames in d.n0.clr with one character values
# gsub is substitute the every occurence of the pattern with the replacement
# pattern is read as: "sa_ followed by one or more other characters (.+)"
# See Regular_expression on wikipedia if you are keen, this is POSIX grep
rownames(d.n0.clr) <- gsub("sa_.+", "S", rownames(d.n0.clr))
rownames(d.n0.clr) <- gsub("td_.+", "T", rownames(d.n0.clr))

# make sure you are using the correct taxon list if you reduce the dataset by cutoff
# as it is now, it is using the entire dataset
# see class code for example
colnames(d.n0.clr) <- gsub(".+__", "", tax)
#palette=palette(c("darkcyan","coral3"))
conds <- data.frame(c(rep(1,length(grep("T", rownames(d.n0.clr)))), rep(2, length(grep("S", rownames(d.n0.clr))))))
colnames(conds) <- "cond"
####

# generate the PCA object
pcx <- prcomp(d.n0.clr)
mvar.clr <- mvar(d.n0.clr)

# plot it
# you control the color of the taxa with the rgb(red,green,blue,transparency) function
# you control the size of the sample and tax labels with the cex() function
coloredBiplot(pcx,col=rgb(0,0,0,0.3),cex=c(0.8,0.4), xlabs.col=conds$cond,
    xlab=paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar.clr, 3)),
    ylab=paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar.clr, 3)),
var.axes=F, arrow.len=0.05, scale=0)
# this gives the proportion of variance explained in the first five components
(pcx$sdev^2/mvar(d.n0.clr))[1:5]

mydata <- as.data.frame(d.n0.clr)
# k-means clustering
# strongly suggests 2 groups
wss <- (nrow(mydata) -1) * sum(apply(mydata,2,var))
# summarize the fit for the first 10 groups
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
        centers=i)$withinss)

layout( matrix(c(1,2,3,4),2,2, byrow=T), widths=c(10,10), height=c(10,10))

par(mar=c(5,4,0,5)+0.1)

# this is the sum of squares distance plot
plot(1:15, wss[1:15], type="b", xlab="Number of Clusters",
    ylab="Win Grp SS")

par(mar=c(2,0,2,0)+0.1)

fit <- kmeans(mydata,2) # fit to 2 clusters
mydata$fit.cluster <- NULL # clear any previous data
mydata <- data.frame(mydata, fit$cluster) # add the cluster data

coloredBiplot(pcx,col=rgb(0,0,0,0.2),cex=c(0.8,0.2),
 xlab=paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar.clr, 3), sep=""),
 ylab=paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar.clr, 3), sep=""),
 xlabs.col=mydata$fit.cluster, arrow.len=0.05, var.axes=F, expand=0.8,  scale=0)

fit <- kmeans(mydata,3) # fit to 2 clusters
mydata$fit.cluster <- NULL # clear any previous data
mydata <- data.frame(mydata, fit$cluster) # add the cluster data

coloredBiplot(pcx,col=rgb(0,0,0,0.2),cex=c(0.8,0.5),
 xlab=paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar.clr, 3), sep=""),
 ylab=paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar.clr, 3), sep=""),
 xlabs.col=mydata$fit.cluster, arrow.len=0.05, var.axes=F, expand=0.8,  scale=0)

fit <- kmeans(mydata,4) # fit to 2 clusters
mydata$fit.cluster <- NULL # clear any previous data
mydata <- data.frame(mydata, fit$cluster) # add the cluster data

coloredBiplot(pcx,col=rgb(0,0,0,0.2),cex=c(0.8,0.2),
 xlab=paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar.clr, 3), sep=""),
 ylab=paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar.clr, 3), sep=""),
 xlabs.col=mydata$fit.cluster, arrow.len=0.05, var.axes=F, expand=0.8,  scale=0)
