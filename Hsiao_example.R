##########
# using the supplementary table from Hsio et al, 2013, Cell Microbiota Modulate Behavioral and Physiological Abnormalities Associated with Neurodevelopmental Disorders
# this shows how I read the table and formatted it for ALDEx2
# d <- read.table("Hsiao_cell/cell7251mmc2/otu_table_S_P_P+Bf.csv", header=T, row.names=1, sep=",",skip=1, comment.char="")
# 
# #put taxonomy into a new variable and then delete
# tax <- d$taxonomy
# d$taxonomy <- NULL
# d.sub <- d[,c.d[1:20]]
# 
# # reduce to OTUs that contain at least 5 counts
# # this reduces to 703 OTUs
# tax.gt5 <- data.frame(tax[which(apply(d.sub,1,sum) > 4)])
# rownames(tax.gt5) <- rownames(d.subgt5)
# d.subgt5 <- data.frame(d.sub[which(apply(d.sub,1,sum) > 4),])
#
# this is a subset of the Hsiao supplementary table dealing only with
# the B. fragilis and its control, this is the table we will use
# in both cases the OTU is the row name
# write.table(d.subgt5, file="hsiao5.txt", sep="\t", quote=F, col.names=NA)
# write.table(tax.gt5, file="tax.txt", sep="\t", quote=F, col.names=NA)
# ##########

# to read the OTU table
d <- read.table("hsiao5.txt", header=T, row.names=1)
tax.d <- read.table("tax.txt", row.names=1, header=T, sep="\t")

library(ALDEx2)

# set up the conditions for comparison
# see the documentation for glm when more than 2 conditions
conds <- c(rep("Bf", 10), rep("C", 10))

# this is the meat of the method
# we generate a set of 128 random Dirichlet instances where the relative
# abundances of each OTU are consistent with the count table and subsequent random
# sampling. Each instance is then clr transformed to place the data in the correct
# space and shape for subsequent analyses

# in general 128 samples are sufficient, and the prior of 0.5 is fine
# a dynamic prior to be used when there are widely divergent read counts per
# sample will be introduced in the next version

x <- aldex.clr(d, mc.samples=128)

# conduct the statistical tests and calculate FDR corrected values
# data are medians of all Dir instances for each OTU
x.t <- aldex.ttest(x, conds)

# calculate the effect sizes for plotting
x.e <- aldex.effect(x, conds)

# merge into one data frame for plotting and examination
x.all <- data.frame(x.t, x.e)

# explore the dataset
aldex.plot(x.all)

diff <- which(x.all$wi.eBH < 0.1)
p.diff <- which(x.all$wi.ep < 0.05)
# or
plot(x.all$diff.win, x.all$diff.btw,pch=19, cex=0.5, xlab="Diff btw", ylab="Diff.win")
# see if anything is below the FDR cutoff
points(x.all$diff.win[diff], x.all$diff.btw[diff],pch=19, cex=0.5, col="red")
# see if anything is below the P cutoff
points(x.all$diff.win[p.diff], x.all$diff.btw[p.diff],pch=19, cex=0.5, col="orange")
# find a particular OTU
points(x.all["145","diff.win"], x.all["145","diff.btw"], cex=1, lwd=2,col="cyan")
# plot the effect size of 1 lines
abline(0,1, col="gray", lty=2)
abline(0,-1, col="gray", lty=2)

#### 
# conclusions
# there are several OTUs that pass a non-FDR threshold of 0.05, however, none
# of them pass when false discovery rate correction is applied
####
# which pass raw P
tax.d[which(x.all$wi.ep < 0.05),]
rownames(d)[which(x.all$wi.ep < 0.05)]

# Example output for a single OTU
# all values are the expected values of the 128 Dirichlet instances
# in general, with sample sizes greater than about 7-8 the Wilcoxon statistic
# will be just as powerful and more reliable than the Welch's test

# x.all["837",]
#         we.ep   we.eBH     wi.ep    wi.eBH   rab.all rab.win.Bf rab.win.C
# 837 0.2413492 0.868727 0.2079292 0.8175367 0.1999454 -0.7850205   1.13357
#     diff.btw diff.win    effect   overlap
# 837 1.783828 3.631399 0.4353015 0.2843751
# 
# we.ep - expected P value from a Welch's t-test  
# we.eBH - expected Benjamini-Hochberg corrected P value
# wi.ep, wi.eBH - same but using Wilcoxon rank test
# 
# rab.all, rab.win.Bf, rab.win.C
# log-ratio abundances relative to the geometric median abundance across all samples, 
# within the Bf samples and within the C samples (see conds). These are log2 based ratios
# so a value of -1 is 2-fold less abundant than the median, a value of 2 is 4-fold more
# abundant than the median, and a value of 0 is equal to the median abundance 
# 
# diff.btw - the difference between the Bf and C values obtained by vector subtraction
# diff.win - the maximum mean absolute difference within the Bf or C group
# 
# effect - the ratio of the btw/win differences
# 
# overlap - the proportion of the time that the Bf and C distributions overlap in the 128 instances
# 
