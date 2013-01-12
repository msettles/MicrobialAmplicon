#### Parse PATRIC files for 16S genes

### Patric Download directory
#seqDir <- "Speciation/Bifidobacteriaceae_patric"
seqDir <- "Speciation/Streptococcaceae_patric"
#Name <- "Bifidobacteriaceae"
Name <- "Streptococcaceae"
align.db <- "Speciation/silva.bacteria.fasta"
source("R/getGenBank.R")

library(Biostrings)
#######
seqFiles <- dir(path=seqDir,pattern=".fna$",full.names=FALSE)
featureFiles <- dir(path=seqDir,pattern=".features.tab$",full.names=FALSE)
seqFiles <- seqFiles[match(sub(".PATRIC.features.tab","",featureFiles),sub(".fna","",seqFiles))]
isolateName <- sub(".PATRIC.features.tab","",featureFiles)
genus <- sapply(strsplit(isolateName,split="_"),"[[", 1L)
species <- sapply(strsplit(isolateName,split="_"),"[[", 2L)
taxonomy <- cbind(genus,species,isolateName)

sequences <- lapply(file.path(seqDir,seqFiles),readDNAStringSet)
annot <- lapply(file.path(seqDir,featureFiles),read.table,sep="\t",header=TRUE,as.is=TRUE,quote="",comment.char="")

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
write.table(anno,file.path(seqDir,paste(Name,"patric.annot.txt",sep=".")),sep="\t",col.names=TRUE,row.names=FALSE)
writeXStringSet(seqs,file.path(seqDir,paste(Name,"patric.fasta",sep=".")))

### CD_HIT REDUNDANCY
system(paste("~/opt/bin/cdhit-est -M 1400 -T 8 -d 200 -c 1.000 -n 9 -i",file.path(seqDir,paste(Name,"patric.fasta",sep=".")),"-o",file.path(seqDir,paste(Name,"patric.reduced.fasta",sep="."))))
#### look at clustering results
cdhit_seq <- readDNAStringSet(file.path(seqDir,paste(Name,"patric.reduced.fasta",sep=".")))
cdhit_cluster <- readLines(file.path(seqDir,paste(Name,"patric.reduced.fasta.clstr",sep=".")))
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

cluster_names <- sapply(lapply(spXcl,function(x) table(as.character(x))),function(x) paste(paste(names(x),x,sep="."),collapse=";"))
cluster_mat_rep <- cluster_mat[is.na(cluster_mat$Identity),]

cluster_mat_rep$sequence_pool[match(names(cluster_names),cluster_mat_rep$Cluster_ID)] <- cluster_names
seqs <- seqs[match(cluster_mat_rep$ID,names(seqs))]

#extraSeqs <- readDNAStringSet("sequence.fasta")

#OUTLIERS <- c(18,19,89,119,164)
#seqo <- seqs[-OUTLIERS]
#cluster_mat_repo <- cluster_mat_rep[-OUTLIERS,]

seqo <- seqs
cluster_mat_repo <- cluster_mat_rep

#writeXStringSet(c(extraSeqs,seqo),"Lactobacillaceae.patric.red.fa")
#write.table(rbind(c(0,width(extraSeqs),names(extraSeqs),NA,NA,0,"Lactobacillus","fornicalis","Lactobacillus crispatus TV1018","fornicalis.1"), cluster_mat_repo),"Lactobacillaceae.patric.red.taxonomy",sep="\t",row.names=F,col.names=T,quote=F)

writeXStringSet(c(seqo),file.path(seqDir,paste(Name,"patric.red.fa",sep=".")))
write.table(cluster_mat_repo,file.path(seqDir,paste(Name,"patric.red.taxonomy",sep=".")),sep="\t",row.names=F,col.names=T,quote=F)
          
########## MOTHUR

system(paste("/Users/mattsettles/opt/bin/mothur \"#align.seqs(candidate=",file.path(seqDir,paste(Name,"patric.red.fa",sep=".")),", template=",align.db,", flip=T, processors=12); filter.seqs(fasta=",file.path(seqDir,paste(Name,"patric.red.align",sep=".")),", processors=12);\"",sep="")) 

Lact.align <- read.table(file.path(seqDir,paste(Name,"patric.red.align.report",sep=".")),sep="\t",header=T,as.is=T)
gb <- get.GenBank(unique(Lact.align$TemplateName))
gb_mapped <- data.frame(ID=names(gb),genus=sapply(strsplit(attr(gb, "species"),split="_"),"[[",1L),species=sapply(strsplit(attr(gb, "species"),split="_"),"[[",2L))

Lact.align <- data.frame(Lact.align,gb_mapped[match(Lact.align$TemplateName,gb_mapped$ID),])

## need to write out this table

#### SEE IF THEY CLUSTER, LOOK FOR OUTLIERS

#### use mothur to produce distance matrix
system(paste("/Users/mattsettles/opt/bin/mothur \"#dist.seqs(fasta=",file.path(seqDir,paste(Name,"patric.red.filter.fasta",sep=".")),", calc=onegap, output=square, processors=12)\"",sep=""))
clusters <- read.table(file.path(seqDir,paste(Name,"patric.red.taxonomy",sep=".")),sep="\t",header=T,as.is=F) 
library(flashClust)
library(WGCNA)
d5k <- read.table(file.path(seqDir,paste(Name,"patric.red.filter.square.dist",sep=".")),skip=1,row.names=1)
hc <- flashClust(as.dist(d5k),method="average")

minModuleSize = 1;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = hc, distM = as.matrix(d5k),method="hybrid",
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

dynamicColors = labels2colors(dynamicMods)
speciesColors = labels2colors(clusters$SPECIES)
genusColors = labels2colors(clusters$GENUS)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file.path(seqDir,paste(Name,"Dendrogram_fullsequence.pdf",sep="_")),width=60,height=16,pointsize=8)
plotDendroAndColors(hc, data.frame(TreeCut=dynamicColors,Species=speciesColors,Genus=genusColors), c("Tree Cut","Species","Genus"),
                    dendroLabels = clusters$GENOME_NAME, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.dendroLabels=1.0,
                    main = paste("Clustering Full length",Name,"16S sequence"))
dev.off()

### Save results
file.copy(from=file.path(seqDir,paste(Name,"patric.red.align",sep=".")),to=file.path("Speciation/SpeciateIT2",paste(Name,"patric.red.align",sep=".")),overwrite=TRUE)
file.copy(from=file.path(seqDir,paste(Name,"patric.red.fa",sep=".")),to=file.path("Speciation/SpeciateIT2",paste(Name,"patric.red.fa",sep=".")),overwrite=TRUE)
file.copy(from=file.path(seqDir,paste(Name,"patric.red.taxonomy",sep=".")),to=file.path("Speciation/SpeciateIT2",paste(Name,"patric.red.taxonomy",sep=".")),overwrite=TRUE)
file.copy(from=file.path(seqDir,paste(Name,"Dendrogram_fullsequence.pdf",sep="_")),to=file.path("Speciation/SpeciateIT2",paste(Name,"Dendrogram_fullsequence.pdf",sep="_")),overwrite=TRUE)

