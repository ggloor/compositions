input <- read.table("td_OTU_tag_mapped_lineage_notech.txt", header=T, row.names=1, check.names=F, comment.char="", sep="\t", skip=1)

library(compositions)
library(ALDEx2)
library(cluster)
library(fpc)
library(MASS)

tax <- input$taxonomy
input$taxonomy <- NULL
t.input<-t(input)

# simplest prior, does not matter much when many OTUS
all <- input+0.5

# samples are column wise
all.clr <- apply(all, 2, function(x){log2(x) - mean(log2(x))})

#build a clr transformed matrix for each individual sample
sampleids_sub1 <- grep("sub1", colnames(all.clr))
sub1.clr <- all.clr[,sampleids_sub1]

sampleids_sub2 <- grep("sub2", colnames(all.clr))
sub2.clr <- all.clr[,sampleids_sub2]

sampleids_sub3 <- grep("sub3", colnames(all.clr))
sub3.clr <- all.clr[,sampleids_sub3]

#make transposed versions of all clr tables
t.all.clr<-t(all.clr)
t.sub1.clr<-t(sub1.clr)
t.sub2.clr<-t(sub2.clr)
t.sub3.clr<-t(sub3.clr)

# get sample-wise euclidian distances for all and individuals
# dist works row-wise requiring the transposed tables
#depending on down stream function need either a matrix, or not
all.clr.dist <- dist(t.all.clr)
sub1.clr.dist<-dist(t.sub1.clr)
sub2.clr.dist<-dist(t.sub2.clr)
sub3.clr.dist<-dist(t.sub3.clr)

m.all.clr.dist <- as.matrix(dist(t.all.clr))
m.sub1.clr.dist<-as.matrix(dist(t.sub1.clr))
m.sub2.clr.dist<-as.matrix(dist(t.sub2.clr))
m.sub3.clr.dist<-as.matrix(dist(t.sub3.clr))
#---------------------------------------------------------------------------
#do kmeans clustering

all.kmean<-kmeans(t.all.clr, 3)

write.table(all.kmean$cluster,"kmean_all_clustered.txt", sep='\t', quote=F, eol='\n')


sub1.kmean<-kmeans(t(sub1.clr),18)
write.table(sub1.kmean$cluster,"kmean_sub1_clustered.txt", sep='\t', quote=F, eol='\n')

sub2.kmean<-kmeans(t(sub1.clr),18)
write.table(sub2.kmean$cluster,"kmean_sub2_clustered.txt", sep='\t', quote=F, eol='\n')

sub3.kmean<-kmeans(t(sub1.clr),18)
write.table(sub3.kmean$cluster,"kmean_sub3_clustered.txt", sep='\t', quote=F, eol='\n')


#plotting from http://stats.stackexchange.com/questions/31083/how-to-produce-a-pretty-plot-of-the-results-of-k-means-cluster-analysis

all.sil<-silhouette(all.kmean$cluster, daisy(t.all.clr)^2)
sub1.sil<-silhouette(sub1.kmean$cluster, daisy(t.sub1.clr)^2)
sub2.sil<-silhouette(sub2.kmean$cluster, daisy(t.sub2.clr)^2)
sub3.sil<-silhouette(sub3.kmean$cluster, daisy(t.sub3.clr)^2)

pdf("kmean_silhouette_plots.pdf")
plot(all.sil)
plot(sub1.sil)
plot(sub2.sil)
plot(sub3.sil)
dev.off()


library(fpc)
pdf("kmean_clusterplot.pdf")
clusplot(t.all.clr, all.kmean$cluster, color=T, shade = T, labels = 0, lines =0)
#clusplot(t(sub1.clr), sub1.kmean$cluster, color=T, shade = T, labels = 0, lines =0)
#clusplot(t(sub2.clr), sub1.kmean$cluster, color=T, shade = T, labels = 0, lines =0)
#clusplot(t(sub3.clr), sub1.kmean$cluster, color=T, shade = T, labels = 0, lines =0)
dev.off()

#find the centers of each of the 3 clusters and find which samples are closest to it
#participant 1 is cluster 3
#participant 2 is cluster 1
#participant 3 is cluster 2
#the output of the kmeans clustering center is a vector of the center OTU clr-transformed values, the next question is which of the samples is closest to it in distance

sub1.center<-all.kmean$centers[3,]
sub2.center<-all.kmean$centers[1,]
sub3.center<-all.kmean$centers[2,]

#push extra column on to the clr tables, rerun the dist function and pull only the distances to the center, sort and find the minimum distance

center.sub1.clr<-as.data.frame(sub1.clr)
center.sub1.clr$center<-sub1.center

center.sub2.clr<-as.data.frame(sub2.clr)
center.sub2.clr$center<-sub2.center

center.sub3.clr<-as.data.frame(sub3.clr)
center.sub3.clr$center<-sub3.center

center.sub1.clr.dist<-as.matrix(dist(t(center.sub1.clr)))
center.sub2.clr.dist<-as.matrix(dist(t(center.sub2.clr)))
center.sub3.clr.dist<-as.matrix(dist(t(center.sub3.clr)))

dist_to_center_sub1<-center.sub1.clr.dist[,"center"]
write.table(sort(dist_to_center_sub1),file="kmeans_sub1_dist_to_center.txt", sep='\t', eol='\n', quote=F)

dist_to_center_sub2<-center.sub2.clr.dist[,"center"]
write.table(sort(dist_to_center_sub2),file="kmeans_sub2_dist_to_center.txt", sep='\t', eol='\n', quote=F)

dist_to_center_sub3<-center.sub3.clr.dist[,"center"]
write.table(sort(dist_to_center_sub3),file="kmeans_sub3_dist_to_center.txt", sep='\t', eol='\n', quote=F)

#for comparisons sake, print out a list of sample distances from the fresh sample

avg_fresh_sub1<-rowMeans(cbind(sub1.clr[,"sub1_NoPres_d0_NAC_ext1_pcr1"], sub1.clr[,"sub1_NoPres_d0_NAC_ext1_pcr2"], sub1.clr[,"sub1_NoPres_d0_NAC_ext2_pcr1"], sub1.clr[,"sub1_NoPres_d0_NAC_ext2_pcr2"])) #calculate the mean CLR values per OTU for the average fresh extracted sample on day 1
fresh.sub1.clr<-as.data.frame(sub1.clr)
fresh.sub1.clr$fresh<-avg_fresh_sub1
fresh.sub1.clr.dist<-as.matrix(dist(t(fresh.sub1.clr)))
dist_to_fresh_sub1<-fresh.sub1.clr.dist[,"fresh"]
write.table(sort(dist_to_fresh_sub1),file="sub1_dist_to_fresh.txt", sep='\t', eol='\n', quote=F)

avg_fresh_sub2<-rowMeans(cbind(sub2.clr[,"sub2_NoPres_d0_NAC_ext1_pcr1"], sub2.clr[,"sub2_NoPres_d0_NAC_ext1_pcr2"], sub2.clr[,"sub2_NoPres_d0_NAC_ext2_pcr1"], sub2.clr[,"sub2_NoPres_d0_NAC_ext2_pcr2"])) #calculate the mean CLR values per OTU for the average fresh extracted sample on day 1
fresh.sub2.clr<-as.data.frame(sub2.clr)
fresh.sub2.clr$fresh<-avg_fresh_sub2
fresh.sub2.clr.dist<-as.matrix(dist(t(fresh.sub2.clr)))
dist_to_fresh_sub2<-fresh.sub2.clr.dist[,"fresh"]
write.table(sort(dist_to_fresh_sub2),file="sub2_dist_to_fresh.txt", sep='\t', eol='\n', quote=F)

avg_fresh_sub3<-rowMeans(cbind(sub3.clr[,"sub3_NoPres_d0_NAC_ext1_pcr1"], sub3.clr[,"sub3_NoPres_d0_NAC_ext1_pcr2"], sub3.clr[,"sub3_NoPres_d0_NAC_ext2_pcr1"], sub3.clr[,"sub3_NoPres_d0_NAC_ext2_pcr2"])) #calculate the mean CLR values per OTU for the average fresh extracted sample on day 1
fresh.sub3.clr<-as.data.frame(sub3.clr)
fresh.sub3.clr$fresh<-avg_fresh_sub3
fresh.sub3.clr.dist<-as.matrix(dist(t(fresh.sub3.clr)))
dist_to_fresh_sub3<-fresh.sub3.clr.dist[,"fresh"]
write.table(sort(dist_to_fresh_sub3),file="sub3_dist_to_fresh.txt", sep='\t', eol='\n', quote=F)




