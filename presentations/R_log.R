library(ALDEx2) # Dirichlet and clr functions
library(compositions) # coloredBiplot function
library(zCompositions) # 0 replacment functions
PATH <- "~/git/compositions/oral/data/"
source("~/git/proprBayes/R/propr-functions.R") # phi association functions
library(igraph) # graph and plotting functions
library(car) # graph and plotting functions

####
# only oral samples were chosen from the entire HMP v35 set
# read in the base dataset for whole mouth
mouth <- read.table(paste(PATH,"mouth_otu.txt", sep=""), header=T,
    row.names=1, sep="\t")

# read in the taxon assignmetns
taxon <- read.table(paste(PATH,"taxon_names.txt", sep=""), header=T,
    row.names=1, sep="\t")

####
# the code that follows is to reduce the dataset to something manageble in
# real time by my laptop. Just about any sensible sample, OTU filtering
# give essentially equivalent answers. Sub-compositional coherence in action!
####

####
# Start biplot plumbing, where data is munged and names changed
####

# remove samples with counts <= 1000
d.col <- mouth[, which(apply(mouth,2,sum) > 1000)]

# remove OTUs with proportion of < 0.0025 in all samples
d.freq <- apply(mouth, 2, function(x){x/sum(x)})
d.bf.0 <- d.col[apply(d.freq, 1, max)>0.0025,]
# align the taxonomy information
tax.0 <- taxon[,1][apply(d.freq, 1, max)>0.0025]

# ensure the taxon is in has a mean abundance of >0.1 counts across the samples
d.bf.0 <- data.frame(d.col[which(apply(d.col, 1, function(x){mean(x) > 0.1})),])
# align the taxonomy information
tax.0 <- taxon[,1][which(apply(d.col, 1, function(x){mean(x) > 0.1}))]

# leaves us with 1446 samples and 5203 OTUs

# replace full names with the single name for display purposes
# grep replaces with the last named taxonomic level
short_tax <- gsub(".+__", "", tax.0)

# make a vector of short ids where the id is a mnemonic for the site
com <- c(rep("T", length(grep("td_.", colnames(d.bf.0))) ),
rep("B", length(grep("bm_", colnames(d.bf.0)))),
rep("A", length(grep("ak_", colnames(d.bf.0)))),
rep("H", length(grep("hp_", colnames(d.bf.0)))),
rep("P", length(grep("pt_", colnames(d.bf.0)))),
rep("S", length(grep("sa_", colnames(d.bf.0)))),
rep("U", length(grep("up_", colnames(d.bf.0)))),
rep("O", length(grep("op_", colnames(d.bf.0)))) )

# make a dataframe of conditions for the coloredBiplot function
conds <- data.frame(c(rep(1,length(grep("td_.", colnames(d.bf.0)))),
    rep(2, length(grep("bm_.", colnames(d.bf.0)))),
    rep(3, length(grep("ak_.", colnames(d.bf.0)))),
    rep(4, length(grep("hp_.", colnames(d.bf.0)))),
    rep(5, length(grep("pt_.", colnames(d.bf.0)))),
    rep(6, length(grep("sa_.", colnames(d.bf.0)))),
    rep(7, length(grep("up_.", colnames(d.bf.0)))),
    rep(8, length(grep("op_.", colnames(d.bf.0))))))
colnames(conds) <- "cond"

# define some colors for each condition
# this could go in the conds dataframe
palette=palette(c("red","darkcyan", "cyan", "darkblue",
    "green", "orange", "magenta", "yellow2"))
####
# End biplot plumbing
####


####
# Start logratio and biplot functions
####

####
# The 0 replacement has an effect on the percentage of variance explained
# CZM: PC1 10.5, PC2 6.1%. Prior: 12.1, 7.0% variance explained
# but there is little if any effect on the overall shape of the projection
# keep in mind that we have thousands of samples and OTUs so either
# is a pretty good result for such a variable dataset
# choose your poison
####

# replace 0 values with an estimate using the count zero multiplicative
# approach from the compositions R package
d.n0 <- cmultRepl(t(d.bf.0), method="CZM", label=0)

# add a prior expectation of 0 to the entire dataset
# naive Bayesian approach
d.n0 <- t(d.bf.0 + 0.5)

# this is a simple clr function
# samples are by row, so in order to preserve this
# t() must wrap the apply function: I hate R
d.n0.clr <- t(apply(d.n0, 1, function(x){log(x) - mean(log(x))}))

# get the metric variance of the dataset, from the compositions package
# this is just the total variance, but it is quicker to compute once and save
mvar.clr <- mvar(d.n0.clr)

# rename the matrix rows and columns (see plumbing)
# warning, this will not work if d.n0.clr is a dataframe!!!
rownames(d.n0.clr) <- com
colnames(d.n0.clr) <- gsub(".+__", "", tax.0)

# this uses the prcomp function from base R
# this is actually kind of magic, because the prcomp function
# makes the PCA from a singular value decomposition -
# which purists argue is the correct way
# the great thing about that is that the resulting PCA is
# defined by a single binary partition
# that is, the clr data is now ILR, so we're legal now
pc.clr <- prcomp(d.n0.clr)

# components in pc.clr$x
# loadings in pc.clr$rotation

# coloredBiplot is from the compositions R package
# most of this is self-explanatory
# xlabs.col is a vector that controls the colors assigned to the
# sample IDs in the order of occurence in the vector and the color palette
# var.axes controls whether arrows are plotted for each OTU
# scale toggles from form (0) to covariance (1) biplot
# form biplots preserve info on samples
# covariance biplots preserve info on OTUs
# to suppress the OTUs from showing: cex=c(0.6,0)
if(print == TRUE) pdf(file="biplot.pdf", height=6, width=6)
coloredBiplot(pc.clr,col=rgb(0,0,0,0.2),cex=c(0.6,0.15),
 xlab=paste("PC1: ", round(sum(pc.clr$sdev[1]^2)/mvar.clr, 3), sep=""),
 ylab=paste("PC1: ", round(sum(pc.clr$sdev[2]^2)/mvar.clr, 3), sep=""),
 xlabs.col=conds$cond, arrow.len=0.05, var.axes=F, expand=0.9,  scale=0)
if(print == TRUE) dev.off()

####
# Here ends the biplot
# for rules of interpretation of compositional biplots
# Aitchison and Greenacre, J. Royal Stat. Society, 2002
####

####
# Start plumbing for tongue vs. cheek comparison
####

# assign tongue and cheek samples to tvsc
tvsc <- data.frame(d.bf.0[,grep("td_", colnames(d.bf.0))],
    d.bf.0[,grep("bm_", colnames(d.bf.0))])

# filter for features that are not 0, and get short taxon names
tvsc.n0 <- tvsc[apply(tvsc, 1, sum) > 0,]
tvsc.tax <-  gsub(".+__", "",tax.0[apply(tvsc, 1, sum) > 0])

conds.a <- c(rep("T", length(grep("td_", colnames(tvsc.n0))) ),
    rep("C", length(grep("bm_", colnames(tvsc.n0)))))
conds.a
####
# End plumbing for tongue vs. cheek
####

# this function is from ALDEx2
# generate 128 Dirichlet Monte Carlo replicates of the dataset
# each replicate contains posterior matrix of probabilities that is consistent
# with the data. We use the observed dataset + 0.5 as the prior
# data must be integers and can include 0
# each MC replicate is clr-transformed
# in essence calculating a distribution of probabilities that are consistent
# with random sampling of the observed dataset
x <- aldex.clr(tvsc.n0, mc.samples=128)

# this function is from ALDEx2
# generate effect size estimates for each of the OTUs
# all values are ratios relative to the geometric mean and not directly abundances
# the results are consistent with the concept of finding those features with the largest
# variance on the simplex that are robust to random sampling of the data
# in practice, features with the greatest error of measurement are close to the low-count
# margin of the simplex
# rab.all: mean ratio abundance
# diff.btw: expected value of the mean signed difference between the two distributions
# diff.win: maximum expected value of the mean unsigned difference within the two distributions
# overlap: proportion of the two distributions that overlaps
# effect: expected value of diff.btw/diff.win
x.e <- aldex.effect(x, conds.a)

# add taxon info to x.e
x.e$tax <- tvsc.tax

# significance tests can be done with aldex.ttest
# this has the same inputs as aldex.effect, and is about 3X slower

# choose effect size bins for display purposes
sig <- which(abs(x.e$effect) >= 1)
sig8 <- which(abs(x.e$effect) >= 0.8 & abs(x.e$effect) < 1)
sig6 <- which(abs(x.e$effect) >= 0.6 & abs(x.e$effect) < .8)
sig4 <- which(abs(x.e$effect) >= 0.4 & abs(x.e$effect) < .6)
sig2 <- which(abs(x.e$effect) >= 0.2 & abs(x.e$effect) < .4)
sig1 <- which(abs(x.e$effect) >= 0.1 & abs(x.e$effect) < .2)

# make the effect plot
if(print == TRUE) pdf("effect.pdf", height=5, width=5)
plot(x.e$diff.win, x.e$diff.btw, pch=19, cex=0.5, col=rgb(0,0,0,0.2), xlab="max group variance", ylab="difference between groups")
points(x.e$diff.win[sig1], x.e$diff.btw[sig1], pch=19, cex=0.6, col=rgb(1,0.4,0,0.5)) # orange
points(x.e$diff.win[sig2], x.e$diff.btw[sig2], pch=19, cex=0.6, col=rgb(0.9,0.7,0,0.5)) # yellowish
points(x.e$diff.win[sig4], x.e$diff.btw[sig4], pch=19, cex=0.6, col=rgb(0,1,0,0.5)) # green
points(x.e$diff.win[sig6], x.e$diff.btw[sig6], pch=19, cex=0.6, col=rgb(0,0.7,1,0.5)) # cyan
points(x.e$diff.win[sig8], x.e$diff.btw[sig8], pch=19, cex=0.6, col=rgb(0,0,1,0.5)) # moderate effect, blue
points(x.e$diff.win[sig], x.e$diff.btw[sig], pch=19, cex=0.6, col=rgb(1,0,0,1)) # largest effect, red
abline(0,-1, lty=2)
abline(0,1, lty=2)
if(print == TRUE) dev.off()

# display the most differential taxa and data
x.e[sig,]

####
# The output from ALDEx2 has a simple interpretation in the simplex and can be visualized
# on a compositional biplot. The effect size provides a measure of the correpondence of the
# variance of the feature to the partition chosen on the simplex. This is explained below.
####

####
# plumbing for the tongue vs cheek biplot
####

# generate the clr matrix and metric variance
clr.tvsx <- apply(tvsc.n0+0.5, 2, function(x){log(x) - mean(log(x))})
mv.tvsx <- mvar(t(clr.tvsx))

# generate the SVD
pcx.tvsx <- prcomp(t(clr.tvsx))

# add letters for body sites
rownames(pcx.tvsx$x) <- conds.a

# make the PCA, it really should be a covariance biplot but the function can't plot the
# variables in color, only the samples
if(print == TRUE) pdf("tc_biplot.pdf", height=5, width=5)
biplot(pcx.tvsx, var.axes=F, cex=c(1,0.01), scale=0,
  col=c("grey", "red"),
  xlab=paste("PC1: ", round(sum(pcx.tvsx$sdev[1]^2)/mv.tvsx,3), sep=""),
  ylab=paste("PC2: ", round(sum(pcx.tvsx$sdev[2]^2)/mv.tvsx,3),sep="")
  )
abline(v=0, lty=2)
abline(h=0, lty=2)
points(pcx.tvsx$rotation[sig1,1],pcx.tvsx$rotation[sig1,2], pch=19, cex=0.5, col=rgb(1,.4,0,1))
points(pcx.tvsx$rotation[sig2,1],pcx.tvsx$rotation[sig2,2], pch=19, cex=0.5, col=rgb(.9,.7,0,1))
points(pcx.tvsx$rotation[sig4,1],pcx.tvsx$rotation[sig4,2], pch=19, cex=0.5, col=rgb(0,1,0,1))
points(pcx.tvsx$rotation[sig6,1],pcx.tvsx$rotation[sig6,2], pch=19, cex=0.5, col=rgb(0,.7,1,1))
points(pcx.tvsx$rotation[sig8,1],pcx.tvsx$rotation[sig8,2], pch=19, cex=0.5, col=rgb(0,0,1,1))
points(pcx.tvsx$rotation[sig,1],pcx.tvsx$rotation[sig,2], pch=19, cex=0.5, col=rgb(1,0,0,1))

)
if(print == TRUE) dev.off()

#####
# Phi as a metric of compositional association
# it combines the correlation of direction of variance, and the corrlation
# of amount of variance into one number. In practice is the standardized ratio between
# two features in a dataset and features are highly correlated when this measure is 0
# see Lovell 2015, PLoS Comp Bio 11:e1004075
#####

# calculate the expected value of phi from Dirichlet Monte-Carlo replicates
# with 0 replacement. This is the same input as was used for the pairwise
# difference by aldex above
# ouput is a dataframe
tvsc.sma.df <- aldex.phi(x)

# choose a cutoff. In practice this can be very low for transcriptomes and higher
# for microbiomes
phi.cutoff <- 0.3

tvsc.sma.lo.phi <- subset(tvsc.sma.df, phi < phi.cutoff)

g <- graph.data.frame(tvsc.sma.lo.phi, directed=FALSE)
g.clust <- clusters(g)
g.df <- data.frame(Systematic.name=V(g)$name, cluster=g.clust$membership,
    cluster.size=g.clust$csize[g.clust$membership])
big <- g.df[which(g.df$cluster.size >= 10),]
colnames(big) <- colnames(g.df)
big <- g.df[which(g.df$cluster.size >= 1),]
colnames(big) <- colnames(g.df)

short_tax <- data.frame(tvsc.tax)
rownames(short_tax) <- rownames(tvsc.n0)

 colours <- c("indianred1", "steelblue3",  "skyblue1", "mediumorchid","royalblue4", "olivedrab3",
   "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0",
    "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF",
   "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC",
   "#669999", "#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00",
   "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED"
)

V(g)$name <- as.vector(short_tax[names(V(g)),])
pdf("phi.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(g)

biplot(pcx.tvsx, var.axes=F, cex=c(1,0.01), scale=0,
  col=c("grey", "red"),
  xlab=paste("PC1: ", round(sum(pcx.tvsx$sdev[1]^2)/mv.tvsx,3), sep=""),
  ylab=paste("PC2: ", round(sum(pcx.tvsx$sdev[2]^2)/mv.tvsx,3),sep="")
  )

#abline(v=0, lty=2)
#abline(h=0, lty=2)

lev <- factor(big$cluster)
for(i in as.numeric(levels(lev))){
nms <- rownames(big)[big$cluster==i]
#print(rownames(big)[big$cluster==i])
#print("")
text(pcx.tvsx$rotation[nms,][,1], pcx.tvsx$rotation[nms,][,2],
    labels = short_tax[rownames(big)[big$cluster==i],1],col=colours[i], cex=0.6)
}
points(pcx.tvsx$rotation[sig,][,1],pcx.tvsx$rotation[sig,][,2], col=rgb(0,0,0,0.5), pch=19, cex=1)

dev.off()
