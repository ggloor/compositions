# read in the entire dataset, this is the metadata
# data files are in v35 download column of: http://hmpdacc.org/HMQCP/
# you want row 10, v35 column ftp
id <- read.table("v35_map_uniquebyPSN.txt", header=TRUE, sep="\t", row.names=1)

# you want row 9, v35 column ftp
# this is the raw otu table
otu <- t( read.table("otu_table_psn_v35.txt", header=T, sep="\t", row.names=1) )
rownames(otu) <- sub("^X", "", rownames(otu))

tax <- otu["Consensus.Lineage",]

sites <- levels(id$HMPbodysubsite)

# grab the data from visit number 1
td <- rownames(id)[which(id$HMPbodysubsite=="Tongue_dorsum" & id$visitno == 1)]
bm <- rownames(id)[which(id$HMPbodysubsite=="Buccal_mucosa" & id$visitno == 1)]
ak <- rownames(id)[which(id$HMPbodysubsite=="Attached_Keratinized_gingiva" & id$visitno == 1)]
hp <- rownames(id)[which(id$HMPbodysubsite=="Hard_palate" & id$visitno == 1)]
pt <- rownames(id)[which(id$HMPbodysubsite=="Palatine_Tonsils" & id$visitno == 1)]
sa <- rownames(id)[which(id$HMPbodysubsite=="Saliva" & id$visitno == 1)]
up <- rownames(id)[which(id$HMPbodysubsite=="Subgingival_plaque" & id$visitno == 1)]
op <- rownames(id)[which(id$HMPbodysubsite=="Supragingival_plaque" & id$visitno == 1)]

# each site goes into its own variable
m1 <- otu[rownames(otu) %in% td,]
m1.n <- apply(m1, 1, function(x){as.numeric(x)})
colnames(m1.n) <- gsub("^700", "td_", rownames(m1))

m2 <- otu[rownames(otu) %in% bm,]
m2.n <- apply(m2, 1, function(x){as.numeric(x)})
colnames(m2.n) <- gsub("^700", "bm_", rownames(m2))

m3 <- otu[rownames(otu) %in% ak,]
m3.n <- apply(m3, 1, function(x){as.numeric(x)})
colnames(m3.n) <- gsub("^700", "ak_", rownames(m3))

m4 <- otu[rownames(otu) %in% hp,]
m4.n <- apply(m4, 1, function(x){as.numeric(x)})
colnames(m4.n) <- gsub("^700", "hp_", rownames(m4))

m5 <- otu[rownames(otu) %in% pt,]
m5.n <- apply(m5, 1, function(x){as.numeric(x)})
colnames(m5.n) <- gsub("^700", "pt_", rownames(m5))

m6 <- otu[rownames(otu) %in% sa,]
m6.n <- apply(m6, 1, function(x){as.numeric(x)})
colnames(m6.n) <- gsub("^700", "sa_", rownames(m6))

m7 <- otu[rownames(otu) %in% up,]
m7.n <- apply(m7, 1, function(x){as.numeric(x)})
colnames(m7.n) <- gsub("^700", "up_", rownames(m7))

m8 <- otu[rownames(otu) %in% op,]
m8.n <- apply(m8, 1, function(x){as.numeric(x)})
colnames(m8.n) <- gsub("^700", "op_", rownames(m8))


#############
# base dataset for whole mouth
mouth <- data.frame(m1.n, m2.n, m3.n, m4.n, m5.n, m6.n, m7.n, m8.n)
taxon <- data.frame(otu["Consensus.Lineage",], colnames(otu))
rownames(taxon) <- colnames(otu)

write.table(mouth, file="~/git/working_papers/oral/mouth_otu.txt", sep="\t", quote=F, col.names=NA)
write.table(taxon, file="~/git/working_papers/oral/taxon_names.txt", sep="\t", quote=F, col.names=NA)
