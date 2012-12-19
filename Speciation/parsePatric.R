#### Parse PATRIC files for 16S genes

seqFiles <- dir(pattern=".fna$")
featureFiles <- dir(pattern=".features.tab$")
seqFiles <- seqFiles[match(sub(".PATRIC.features.tab","",featureFiles),sub(".fna","",seqFiles))]
isolateName <- sub(".PATRIC.features.tab","",featureFiles)
genus <- sapply(strsplit(isolateName,split="_"),"[[", 1L)
species <- sapply(strsplit(isolateName,split="_"),"[[", 2L)
taxonomy <- cbind(genus,species,isolateName)

sequences <- lapply(seqFiles,readDNAStringSet)
annot <- lapply(featureFiles,read.table,sep="\t",header=TRUE,as.is=TRUE,quote="",comment.char="")

has16s <- sapply(annot,function(x) length(intersect(grep("Small Subunit",x$PRODUCT),which(x$FEATURE_TYPE=="rRNA"))))
sequence <- sequences[has16s>0]
annot <- annot[has16s>0]
taxonomy <- taxonomy[has16s>0,]
taxonomy <- cbind(taxonomy,has16s[has16s>0])

annot <- lapply(annot,function(x) x[intersect(grep("Small Subunit",x$PRODUCT),which(x$FEATURE_TYPE=="rRNA")),]) 

yankSequence <- function(seqs,anno){
  tmp <- apply(anno,1,function(x){ seq <- subseq(x=seqs[grep(x["ACCESSION"], names(seqs))],start=as.numeric(x["START"]),end=as.numeric(x["END"])); if(as.vector(x["STRAND"] == "+")){seq <- reverseComplement(seq)};seq})
  DNAStringSet(sapply(tmp,as.character))
}

seqs <- sapply(seq.int(1,length(annot)),function(x) yankSequence(sequence[[x]],annot[[x]]))
seqs <- DNAStringSet(unlist(sapply(seqs,as.character)))

names(seqs) <- unlist(sapply(annot,function(anno) paste(anno$ACCESSION,":",anno$START,"-",anno$END,":",anno$STRAND, ":",gsub(" ","_",anno$GENOME_NAME),sep="")))
anno <- do.call("rbind", annot)
longseqs <- which(width(seqs) > 1450)
seqs <- seqs[longseqs]
anno <- anno[longseqs,]
anno$SEQID <- paste(anno$ACCESSION,":",anno$START,"-",anno$END,":",anno$STRAND, ":",gsub(" ","_",anno$GENOME_NAME),sep="")
anno$GENUS <- sapply(strsplit(anno$GENOME_NAME,split=" "),"[[",1L)
anno$SPECIES <- sapply(strsplit(anno$GENOME_NAME,split=" "),"[[",2L)
write.table(anno,"Lactobacillaceae.patric.annot.txt",sep="\t",col.names=TRUE,row.names=FALSE)
writeXStringSet(seqs,"Lactobacillaceae.patric.fasta")


### CD_HIT REDUNDANCY
system("~/opt/bin/cdhit-est -M 1400 -T 8 -d 200 -c 1.000 -n 9 -i Lactobacillaceae.patric.fasta -o Lactobacillaceae.patric.reduced.fasta")
#### look at clustering results
cdhit_seq <- readDNAStringSet("Lactobacillaceae.patric.reduced.fasta")
cdhit_cluster <- readLines("Lactobacillaceae.patric.reduced.fasta.clstr")
clust <- grep("^>Cluster",c(cdhit_cluster,">Cluster"))
cdhit_cluster <- paste(cdhit_cluster,"C",rep(seq.int(1,length(diff(clust))),times=diff(clust)),sep=" ")
cdhit_cluster <- cdhit_cluster[-clust]
test_split <- strsplit(cdhit_cluster,split="\t|nt, >|\\.\\.\\. at +/|\\.\\.\\. |%C |C ")
cluster_mat <- data.frame(matrix(unlist(test_split),ncol=5,byrow=T),stringsAsFactors=F)
cluster_mat$X4 <- as.numeric(gsub("at [+]\\/|% ","",cluster_mat$X4))
cluster_mat$errors <- as.numeric(cluster_mat$X2) - as.numeric(cluster_mat$X2)*(cluster_mat$X4/100)

cluster_mat[is.na(cluster_mat$errors),"errors"] <- 0
cluster_mat <- cluster_mat[order(as.numeric(cluster_mat$X5),-is.na(cluster_mat$X4)),]
colnames(cluster_mat) <- c("ord","Len","ID","Identity","Cluster_ID","Error")


cluster_mat <- cluster_mat[match(names(seqs),cluster_mat$ID),]

cluster_mat <- data.frame(cluster_mat,anno[,c("GENUS","SPECIES","GENOME_NAME")])

spXcl<- split(cluster_mat$SPECIES, cluster_mat$Cluster_ID)
spXcl.len <- sapply(spXcl,function(x) length(unique(x)))

cluster_names <- sapply(sapply(spXcl,function(x) table(as.character(x))),function(x) paste(paste(names(x),x,sep="."),collapse=";"))
cluster_mat_rep <- cluster_mat[is.na(cluster_mat$Identity),]

cluster_mat_rep$sequence_pool[match(names(cluster_names),cluster_mat_rep$Cluster_ID)] <- cluster_names
seqs <- seqs[match(cluster_mat_rep$ID,names(seqs))]

extraSeqs <- readDNAStringSet("sequence.fasta")

OUTLIERS <- c(18,19,89,119,164)
seqo <- seqs[-OUTLIERS]
cluster_mat_repo <- cluster_mat_rep[-OUTLIERS,]
writeXStringSet(c(extraSeqs,seqo),"Lactobacillaceae.patric.red.fa")

write.table(rbind(c(0,width(extraSeqs),names(extraSeqs),NA,NA,0,"Lactobacillus","fornicalis","Lactobacillus crispatus TV1018","fornicalis.1"), cluster_mat_repo),"Lactobacillaceae.patric.red.taxonomy",sep="\t",row.names=F,col.names=T,quote=F)

########## MOTHUR

system("/Users/mattsettles/opt/bin/mothur \"#align.seqs(candidate=Lactobacillaceae.patric.red.fa, template=../silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=Lactobacillaceae.patric.red.align, processors=12);\"") 

Lact.align <- read.table("Lactobacillaceae.patric.red.align.report",sep="\t",header=T,as.is=T)
source("../../R/getGenBank.R")
gb <- get.GenBank(unique(Lact.align$TemplateName))
gb_mapped <- data.frame(ID=names(gb),genus=sapply(strsplit(attr(gb, "species"),split="_"),"[[",1L),species=sapply(strsplit(attr(gb, "species"),split="_"),"[[",2L))

Lact.align <- data.frame(Lact.align,gb_mapped[match(Lact.align$TemplateName,gb_mapped$ID),])


#### SEE IF THEY CLUSTER, LOOK FOR OUTLIERS

#### use mothur to produce distance matrix
system("/Users/mattsettles/opt/bin/mothur \"#dist.seqs(fasta=Lactobacillaceae.patric.red.filter.fasta, calc=onegap, output=square, processors=12)\"")
clusters <- read.table("Lactobacillaceae.patric.red.taxonomy",sep="\t",header=T,as.is=F) 
library(flashClust)
library(WGCNA)
d5k <- read.table("Lactobacillaceae.patric.red.filter.square.dist",skip=1,row.names=1)
hc <- flashClust(as.dist(d5k),method="single")

minModuleSize = 1;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = hc, distM = as.matrix(d5k),method="hybrid",
                            deepSplit = TRUE, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

dynamicColors = labels2colors(dynamicMods)
speciesColors = labels2colors(clusters$SPECIES)
genusColors = labels2colors(clusters$GENUS)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf("Dendrogram_fullsequence.pdf",width=24,height=8,pointsize=8)
plotDendroAndColors(hc, data.frame(TreeCut=dynamicColors,Species=speciesColors,Genus=genusColors), c("Tree Cut","Species","Genus"),
                    dendroLabels = clusters$SPECIES, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.dendroLabels=0.5,
                    main = "Clustering Full length Lactobacillaceae 16S sequence")
dev.off()
