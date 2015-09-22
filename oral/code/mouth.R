library(ALDEx2)
library(compositions)
library(zCompositions)

#############
# read in the base dataset for whole mouth

mouth <- read.table("~/git/compositions/oral/data/mouth_otu.txt", header=T, row.names=1, sep="\t")
taxon <- read.table("~/git/compositions/oral/data/taxon_names.txt", header=T, row.names=1, sep="\t")

# ids have leading 700 stripped off
# samples named as follows:
# td_ Tongue Dorsum
# bm_ Buccal mucosa
# ak_ Attached Keratinzed Gingiva
# hp_ Hard Palate
# pt_ Palatine Tonsils
# sa_ Saliva
# up_ Subgingival Plaque (under plaque)
# op_ Supragingival Plaque (over plaque)

#######
# ids in same order in taxon and in mouth
# mouth.n0 <- mouth[which(apply(mouth, 1, sum) > 0),]
# all rows have at least one count! shit can't reduce this

# remove samples with < 1000 reads
d.col <- mouth[, which(apply(mouth,2,sum) > 1000)]

# First look at the entire mouth at the OTU level
# reduce the dataset to those that are present in at least x% of the samples
# have 1101 of 45383 OTUs left at 20%
# samples are now in row, OTUs in columns

# cutoff is X percent
cutoff = 1-0.50
d.bf.0 <- data.frame(d.col[which(apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) < cutoff),])
tax.0 <- taxon[,1][which(apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) < cutoff)]

# replace 0 values with an estimated value (zCompositions package) this is slow
d.bf <-cmultRepl(t(d.bf.0),  label=0)

# g full names with the last taxonomic name
short_tax <- gsub(".+__", "", tax.0)

# replace sample names with short ids
com <- c(rep("T", length(grep("td_.", colnames(d.bf.0))) ),
rep("B", length(grep("bm_", colnames(d.bf.0)))),
rep("A", length(grep("ak_", colnames(d.bf.0)))),
rep("H", length(grep("hp_", colnames(d.bf.0)))),
rep("P", length(grep("pt_", colnames(d.bf.0)))),
rep("S", length(grep("sa_", colnames(d.bf.0)))),
rep("U", length(grep("up_", colnames(d.bf.0)))),
rep("O", length(grep("op_", colnames(d.bf.0)))) )

##########
# this uses the princomp function as per the compositions package instructions

# acomp expects samples in rows
bi <- acomp(d.bf)
# No. corrected values:  944574
# so we have to be careful interpreting this, we probably want to filter more aggressively

rownames(bi) <- com
names(bi) <- short_tax

pcx <- princomp(bi)

conds <- data.frame(c(rep(1,length(grep("td_.", colnames(d.bf.0)))), rep(2, length(grep("bm_.", colnames(d.bf.0)))), rep(3, length(grep("ak_.", colnames(d.bf.0)))), rep(4, length(grep("hp_.", colnames(d.bf.0)))), rep(5, length(grep("pt_.", colnames(d.bf.0)))), rep(6, length(grep("sa_.", colnames(d.bf.0)))), rep(7, length(grep("up_.", colnames(d.bf.0)))), rep(8, length(grep("up_.", colnames(d.bf.0))))))
colnames(conds) <- "cond"

palette=palette(c("red","cyan", "darkcyan", "darkblue", "green", "orange", "magenta", "yellow2"))
coloredBiplot(pcx,col="black",cex=c(0.6,0.5), xlabs.col=conds$cond, arrow.len=0.05, var.axes=F, expand=0.9,  scale=0)
sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)

#############

#############
# this uses the prcomp function that can be used for high-dimensional data ala MASS
# in the mouth data, they are very similar

dev.new()
d.n0 <- cmultRepl(t(d.bf.0), label=0, method="CZM")

d.n0.clr <- apply(d.n0, 2, function(x){log(x) - mean(log(x))})

rownames(d.n0.clr) <- com
colnames(d.n0.clr) <- gsub(".+__", "", tax.0)
pc.clr <- prcomp(d.n0.clr)
biplot(pc.clr, cex=c(0.5,0.5), scale=0, arrow.len=0, var.axes=F)
coloredBiplot(pc.clr,col="black",cex=c(0.6,0.5), xlabs.col=conds$cond, arrow.len=0.05, var.axes=F, expand=0.9,  scale=0)
sum(pc.clr$sdev[1]^2)/mvar(d.n0.clr)
sum(pc.clr$sdev[2]^2)/mvar(d.n0.clr)
sum(pc.clr$sdev[3]^2)/mvar(d.n0.clr)

##############

sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)

########################
# Tongue vs Palatine Tonsils mucosa : hard to see a difference 18% explained
# Tongue vs Buccal mucosa: clear separation 22% variance
# Tongue vs Hard Palate: poor separation 16% variance
# Tongue vs Saliva: poor separation 17% variance
# Tongue vs Attached gingiva: clear separation 30% variance
# Tongue vs supra plaque: clear separation 29% variance
# Tongue vs sub plaque: clear separation 27% variance
# otu level

TB <- data.frame(mouth[,grep("up_", colnames(mouth))], mouth[,grep("op_", colnames(mouth))])

# remove samples with < 1000 reads
d.col <- TB[, which(apply(TB,2,sum) > 1000)]

# reduce the dataset to those that are present in at least 40% of the samples
# have 1101 of 45383 OTUs left at 20%
# samples are now in row, OTUs in columns

# cutoff is X percent 40% not 0 is the max
cutoff = 1-0.50
d.bf.0 <- data.frame(d.col[which(apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) < cutoff),])
tax.0 <- taxon[,1][which(apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) < cutoff)]

# replace 0 values with an estimated value (zCompositions package) this is slow
d.bf <-cmultRepl(t(d.bf.0),  label=0)

# acomp expects samples in rows
bi <- acomp(d.bf)
# No. corrected values:  944574
# so we have to be careful interpreting this, we probably want to filter more aggressively

#replace full names with the last taxonomic name
names(bi) <- gsub(".+__", "", tax.0)

com <- c(rep("T", length(grep("up_.", colnames(d.bf.0))) ),
rep("S", length(grep("op_", colnames(d.bf.0))))
 )
rownames(bi) <- com

pcx <- princomp(bi)

# this can be very slow with all the grepping
conds <- data.frame(c(rep(1,length(grep("up_.", colnames(d.bf.0)))), rep(2, length(grep("op_.", colnames(d.bf.0)))) )) #, rep(3, length(grep("op_.", colnames(d.bf.0)))), rep(4, length(grep("hp_.", colnames(d.bf.0)))), rep(5, length(grep("pt_.", colnames(d.bf.0)))), rep(6, length(grep("op_.", colnames(d.bf.0)))), rep(7, length(grep("op_.", colnames(d.bf.0)))), rep(8, length(grep("op_.", colnames(d.bf.0))))))
colnames(conds) <- "cond"

palette=palette(c("red","cyan", "darkcyan", "darkblue", "green", "orange", "magenta", "yellow2"))
# scale = 0 is a form biplot - you scale by the arrows
# scale = 1 is a covariance biplot - you scale by the samples
coloredBiplot(pcx,col="black",cex=c(0.6,0.5), xlabs.col=conds$cond, arrow.len=0.05, var.axes=T, expand=0.9,  scale=0)

sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)

########################
# Then look at the entire mouth at the lowest named level (mainly genus and family)
mouth.agg <- aggregate(mouth, by=list(taxon[,1]), FUN=sum)
tax.agg <- mouth.agg$Group.1
mouth.agg$Group.1 <- NULL

d.col <- mouth.agg[, which(apply(mouth.agg,2,sum) > 1000)]
rownames(d.col) <- tax.agg

d.bf.0 <- data.frame(d.col[apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) <0.90,])

d.bf <-cmultRepl(t(d.bf.0),  label=0)
# acomp expects samples in rows
bi <- acomp(d.bf)

#replace full names with the last taxonomic name
names(bi) <- gsub(".+__", "", names(bi))

com <- c(rep("T", length(grep("td_.", colnames(d.bf.0))) ),
rep("B", length(grep("bm_", colnames(d.bf.0)))),
rep("A", length(grep("ak_", colnames(d.bf.0)))),
rep("H", length(grep("hp_", colnames(d.bf.0)))),
rep("P", length(grep("pt_", colnames(d.bf.0)))),
rep("S", length(grep("sa_", colnames(d.bf.0)))),
rep("U", length(grep("up_", colnames(d.bf.0)))),
rep("O", length(grep("op_", colnames(d.bf.0)))) )
rownames(bi) <- com
pcx <- princomp(bi)

conds <- data.frame(c(rep(1,length(grep("td_.", colnames(d.bf.0)))), rep(2, length(grep("bm_.", colnames(d.bf.0)))), rep(3, length(grep("ak_.", colnames(d.bf.0)))), rep(4, length(grep("hp_.", colnames(d.bf.0)))), rep(5, length(grep("pt_.", colnames(d.bf.0)))), rep(6, length(grep("sa_.", colnames(d.bf.0)))), rep(7, length(grep("up_.", colnames(d.bf.0)))), rep(8, length(grep("up_.", colnames(d.bf.0))))))
colnames(conds) <- "cond"

palette=palette(c("darkcyan","cyan", "red", "darkblue", "green", "orange", "magenta", "yellow2"))
coloredBiplot(pcx,col="black",cex=c(0.6,0.5), xlabs.col=conds$cond, arrow.len=0.05, var.axes=T, expand=0.9,  scale=0)

dim(d.bf.0)
sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)


#############
# plaque only
plaque <- data.frame(m7.n, m8.n)
plaque.agg <- aggregate(plaque, by=list(taxon[,1]), FUN=sum)
tax.agg <- plaque.agg$Group.1
plaque.agg$Group.1 <- NULL

d.col <- plaque.agg[, which(apply(plaque.agg,2,sum) > 1000)]
rownames(d.col) <- tax.agg

d.bf.0 <- data.frame(d.col[apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) <0.70,])
# 68 taxa in 619 samples


d.bf <-cmultRepl(t(d.bf.0),  label=0)
# acomp expects samples in rows
bi <- acomp(d.bf)

#replace full names with the last taxonomic name
names(bi) <- gsub(".+__", "", names(bi))


com <- c(rep("U", length(grep("up_", colnames(d.bf.0)))),
rep("O", length(grep("op_", colnames(d.bf.0)))) )

rownames(bi) <- com
pcx <- princomp(bi)

conds <- data.frame(c(rep(7, length(grep("up_.", colnames(d.bf.0)))), rep(8, length(grep("op_.", colnames(d.bf.0))))))
colnames(conds) <- "cond"

palette=palette(c("magenta", "black"))
coloredBiplot(pcx,col="black",cex=0.4, xlabs.col=conds$cond, arrow.len=0.05)

dim(d.bf.0)
sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)

con <- c(rep(1, length(grep("up_.", colnames(d.bf.0)))), rep(2, length(grep("op_.", colnames(d.bf.0)))))

x <- aldex.clr(d.bf.0)
x.e <- aldex.effect(x,con)
x.t <- aldex.ttest(x,con)
x.all <- data.frame(x.e,x.t)

#############
# Attached gingiva vs saliva (AS)
AS <- data.frame(m3.n, m6.n)

# ids in same order in taxon and in AP
# all rows have at least one count! shit
# AS.n0 <- AP[which(apply(AP, 1, sum) > 0),]

AS.agg <- aggregate(AS, by=list(taxon[,1]), FUN=sum)
tax.agg <- AS.agg$Group.1
AS.agg$Group.1 <- NULL

d.col <- AS.agg[, which(apply(AS.agg,2,sum) > 1000)]
rownames(d.col) <- tax.agg

d.bf.0 <- data.frame(d.col[apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) <0.90,])

d.bf <-cmultRepl(t(d.bf.0),  label=0)
# acomp expects samples in rows
bi <- acomp(d.bf)

#replace full names with the last taxonomic name
names(bi) <- gsub(".+__", "", names(bi))

com <- c(rep("A", length(grep("ak_", colnames(d.bf.0)))),
rep("S", length(grep("sa_", colnames(d.bf.0)))) )
rownames(bi) <- com
pcx <- princomp(bi)

conds <- data.frame(c(rep(3, length(grep("ak_.", colnames(d.bf.0)))), rep(5, length(grep("sa_.", colnames(d.bf.0)))))
colnames(conds) <- "cond"

palette=palette(c("red","orange"))
coloredBiplot(pcx,col="black",cex=c(0.6,0.5), xlabs.col=conds$cond, arrow.len=0.05, var.axes=T, expand=1.2,  scale=0)

dim(d.bf.0)
sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)
library(ALDEx2)
con <- c(rep(1, length(grep("ak_.", colnames(d.bf.0)))), rep(2, length(grep("sa_.", colnames(d.bf.0)))))

x <- aldex.clr(d.bf.0)
x.e <- aldex.effect(x,con)
x.t <- aldex.ttest(x,con)
x.all <- data.frame(x.e,x.t)

#############
# Attached tongue vs Palate (TH)
TH <- data.frame(m1.n, m4.n)

# ids in same order in taxon and in AP
# all rows have at leTHt one count! shit
# TH.n0 <- AP[which(apply(AP, 1, sum) > 0),]

TH.agg <- aggregate(TH, by=list(taxon[,1]), FUN=sum)
tax.agg <- TH.agg$Group.1
TH.agg$Group.1 <- NULL

d.col <- TH.agg[, which(apply(TH.agg,2,sum) > 1000)]
rownames(d.col) <- tax.agg

d.bf.0 <- data.frame(d.col[apply(d.col, 1, function(x){length(which(x == 0))/length(x)}) <0.90,])

# 68 taxa in 619 samples


d.bf <-cmultRepl(t(d.bf.0),  label=0)
# acomp expects samples in rows
bi <- acomp(d.bf)

#replace full names with the lTHt taxonomic name
names(bi) <- gsub(".+__", "", names(bi))


com <- c(rep("T", length(grep("td_", colnames(d.bf.0)))),
rep("H", length(grep("hp_", colnames(d.bf.0)))) )
rownames(bi) <- com
pcx <- princomp(bi)

conds <- data.frame(c(rep(1, length(grep("td_.", colnames(d.bf.0)))), rep(2, length(grep("hp_.", colnames(d.bf.0))))))
colnames(conds) <- "cond"

palette=palette(c("darkcyan","red"))
coloredBiplot(pcx,col="black",cex=c(0.6,0.5), xlabs.col=conds$cond, arrow.len=0.05, var.axes=T, expand=1.2,  scale=0)

dim(d.bf.0)
sum(pcx$sdev[1]^2)/mvar(bi)
sum(pcx$sdev[2]^2)/mvar(bi)
sum(pcx$sdev[3]^2)/mvar(bi)

con <- c(rep(1, length(grep("td_.", colnames(d.bf.0)))), rep(2, length(grep("hp_.", colnames(d.bf.0)))))

x <- aldex.clr(d.bf.0)
x.e <- aldex.effect(x,con)
x.t <- aldex.ttest(x,con)
x.all <- data.frame(x.e,x.t)
