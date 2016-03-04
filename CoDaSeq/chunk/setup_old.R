## ---- setup ----

# we need these libraries for the colored biplot, and for the 0 replacement function
library(compositions)
library(zCompositions)
source("chunk/codaMB_functions.R")
# read the table, with column and row names, columns tab delimited
# samples are by column, variables are by row
d.bf.1 <- read.table("tongue_vs_cheek.txt", header=T, row.names=1, sep="\t")

# move taxonomy info to a vector
tax.0 <- d.bf.1$tax
# remove the taxonomy column
d.bf.1$tax <- NULL

# keep only those samples with > min.reads
# the which command returns a vector of numbers corresponding to the columns that have
# > n reads. Try it on its own in the console
min.reads <- 5000
 d.bf.0 <- d.bf.1[,which(apply(d.bf.1,2,sum) > min.reads)]

# simplify the dataset, by keeping only taxa present in x fraction of the samples
# the which function returns a list of row numbers where the number of 0 values in the
# row is greater than the cutoff fraction
cutoff = .3
 d.subset <- data.frame(d.bf.0[which(apply(d.bf.0, 1,
    function(x){length(which(x != 0))/length(x)}) > cutoff),])
# d.subset <- codamb.filter(d.bf.1, min.reads=min.reads, min.fraction=cutoff, min.prop=0,
#    samples.by.row=FALSE)
#tax.subset <- tax.0[which(apply(d.bf.0, 1,
#    function(x){length(which(x != 0))/length(x)}) > cutoff)]
# oh boy, so much here to explain
# basically I am generating a shortened name for each OTU for display purposes
# you could do this by hand if you don't know grep
# two ways to substitute row and column labels
# gsub replaces any number of characters denoted by the "." that are followed by
# two underscores with nothing from the list of taxa names in the reduced dataset

# short_tax <- gsub(".+__", "", tax.subset)

# com is being assigned a vector of values where the "T" is replicated the number of times
# that we observe td_ in the column names of the reduced data subset
# com <- c(rep("T", length(grep("td_.", colnames(d.subset))) ),
#    rep("B", length(grep("bm_.", colnames(d.subset))) ))

