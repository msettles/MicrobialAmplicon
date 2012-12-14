
## Update metadata
### Retreive and Send metadata or mapping data
library(RSQLite)
library(Biostrings)
microbe.amplicon.home <- "/mnt/home/msettles/CodeProjects/Rpackages/MicrobialAmplicon"
source(file.path(microbe.amplicon.home,"R","DButilities.R"))


basedir <- "/mnt/lfs/msettles/projects/Amplicon_Preprocessing"

con <- dbCon(file.path(basedir,"amplicondataV2.0.sqlite"))

#######################################################################################
## Speciation
#######################################################################################
#con <- dbCon()
## Lactobacillus
## extract genus reads and write to file
project <- "Witkin-VVS"
output.genus.reads(con,project=project,genus="Lactobacillus")

nproc=8
mothur.template="/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta" 

## cluster using cdhit require 99.5% identity
system("cdhit-est -T 8 -d 0 -c 0.995 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced.0.995.fasta")
# 
## look at clustering results
cdhit_cluster <- readLines("454Reads.Lactobacillus.reduced.0.995.fasta.clstr")
clust <- grep("^>Cluster",c(cdhit_cluster,">Cluster"))
cdhit_cluster <- paste(cdhit_cluster,"C",rep(seq.int(1,length(diff(clust))),times=diff(clust)),sep=" ")
cdhit_cluster <- cdhit_cluster[-clust]
test_split <- strsplit(cdhit_cluster,split="\t|nt, >|\\.\\.\\. at +/|\\.\\.\\. |%C |C ")
cluster_mat <- data.frame(matrix(unlist(test_split),ncol=5,byrow=T),stringsAsFactors=F)
cluster_mat$X4 <- as.numeric(gsub("at [+]\\/|% ","",cluster_mat$X4))
cluster_mat$errors <- as.numeric(cluster_mat$X2) - as.numeric(cluster_mat$X2)*(cluster_mat$X4/100)

cluster_mat[is.na(cluster_mat$errors),"errors"] <- 0
cluster_mat <- cluster_mat[order(as.numeric(cluster_mat$X5),-is.na(cluster_mat$X4)),]

system("mothur \"#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.995.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.995.align, processors=12);\"")
system("mv 454Reads.filter 454Reads.0.995.filter")

#### TEST WHICH c value (cluster percent identity is best)
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 1.000 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-1.0.fasta &
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 0.995 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-0.995.fasta &
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 0.990 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-0.990.fasta &
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 0.985 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-0.985.fasta &
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 0.980 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-0.980.fasta &
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 0.975 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-0.975.fasta &
# nohup cdhit-est -M 1400 -T 8 -d 0 -c 0.970 -n 9 -i 454Reads.Lactobacillus.fasta -o 454Reads.Lactobacillus.reduced-0.970.fasta &
# 
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.1.0.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.1.0.align, processors=12);"
# mv 454Reads.filter 454Reads.1.0.filter
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.995.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.995.align, processors=12);"
# mv 454Reads.filter 454Reads.0.995.filter
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.990.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.990.align, processors=12);" 
# mv 454Reads.filter 454Reads.0.990.filter
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.985.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.985.align, processors=12);" 
# mv 454Reads.filter 454Reads.0.985.filter
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.980.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.980.align, processors=12);" 
# mv 454Reads.filter 454Reads.0.980.filter
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.975.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.975.align, processors=12);" 
# mv 454Reads.filter 454Reads.0.975.filter
# mothur "#align.seqs(candidate=454Reads.Lactobacillus.reduced.0.970.fasta, template=/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta, flip=T, processors=12); filter.seqs(fasta=454Reads.Lactobacillus.reduced.0.970.align, processors=12);" 
# mv 454Reads.filter 454Reads.0.970.filter

## names, no longer needed after rewrite
## names <- cbind(cluster_mat$X3[!duplicated(cluster_mat$X5)],tapply(cluster_mat$X3,cluster_mat$X5,paste,collapse=","))
## write.table(names,file="454Reads.Lactobacillus.reduced.fasta.clstr2",sep=",",row.names=F,col.names=F,quote=F)

## add in Training species sequences
## maybe read them into R and then write out
# speciateIThome <- "/mnt/home/msettles/opt/speciateit"
# speciesFile <- file.path(speciateIThome,"spp-data/Lactobacillus/rdp_vagi1200Len_Lactobacillus.V3V1.550bp_nr.fa")
# system(paste("cat ",speciesFile, " 454Reads.Lactobacillus.reduced.fasta > 454Reads.Lactobacillus.reduced.combined.fa",sep=""))

## Prepare the reference

filter <- readLines("454Reads.0.995.filter")

ref_seqs <- readBStringSet("Lactobacillaceae.patric.red.align")
keep <- which(unlist(strsplit(filter,"")[[1]]==1))
trimmed_ref <- DNAStringSet(gsub("[.]|-","",sapply(strsplit(as.character(ref_seqs),""),function(x) paste( x[keep] , collapse=""))))

#names(trimmed_ref) <- sapply(strsplit(names(trimmed_ref),split="-"),"[[",1L)
names(trimmed_ref) <-  paste("S",seq.int(1,length(trimmed_ref)),sep="")
writeXStringSet(c(trimmed_ref,readDNAStringSet("454Reads.Lactobacillus.reduced.0.995.fasta")),"454Reads.Lactobacillus.reduced.combined.fasta")

spp_taxon <- read.table("Lactobacillaceae.patric.red.taxonomy",sep="\t",header=TRUE,as.is=TRUE)
spp_taxon <- data.frame(ID=paste("S",seq.int(1,length(trimmed_ref)),sep=""),Species=paste(spp_taxon[,"GENUS"],spp_taxon[,"SPECIES"],sep=" "))
write.table(spp_taxon,file="Lactobacillaceae.patric.taxon",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
#spp_taxon <- file.path(speciateIThome,"spp-data/Lactobacillus/rdp_vagi1200Len_Lactobacillus.V3V1.550bp_nr.taxon")

## Use mothur to compute the Distance between sequences 
#### Align the seuqences to silva
#### Filter the alignment
#### Compute distance matrix

system(paste("mothur \"#align.seqs(candidate=454Reads.Lactobacillus.reduced.combined.fasta, template=", mothur.template ,", flip=T, processors=",nproc,"); ",
              "filter.seqs(fasta=454Reads.Lactobacillus.reduced.combined.align, processors=",nproc,"); ",
              "dist.seqs(fasta=454Reads.Lactobacillus.reduced.combined.filter.fasta, calc=onegap, output=square, processors=",nproc,");\"",sep=""))

#system(paste("mothur \"#dist.seqs(fasta=454Reads.Lactobacillus.reduced.0.995.filter.fasta, calc=onegap, output=square, processors=",nproc,");\"",sep=""))

dist.mat <- "454Reads.Lactobacillus.reduced.combined.filter.square.dist"
outFile <- "454Reads.Lactobacillus.reduced.combined.filter.sing.hclust.membStr"

### Use FlashClust and Vicut to produce clusters
library(flashClust)
d5k <- read.table(dist.mat,skip=1,row.names=1)
hc <- flashClust(as.dist(d5k),method="single")

# using rowIds for leaves
### produces file ready for vicut
hclust2merges2 <- function(hc,rowIds,filename)
{
  sink(file=filename)
  internal_id = -1
  for (i in 1:nrow(hc$merge))
    {
      # Less than 1 implies a singleton. Reverse all the numbers.
      u = hc$merge[i,1]*(-1)
      v = hc$merge[i,2]*(-1)

      if ( u > 0 ) u <- rowIds[u]
      if ( v > 0 ) v <- rowIds[v]

      cat(u,v,internal_id,"\n",file=filename,sep="\t",append=TRUE)
      internal_id = internal_id - 1
    }
  sink()
}

hclust2merges2(hc,rownames(d5k),outFile)


system(paste("vicut -t ", outFile, "  -a Lactobacillaceae.patric.taxon -o 454Reads.Lactobacillus.reduced.combined.filter.sing.hclust.minNodeCut.dir",sep=""))

### Read in clusters computed by vicut, annotate by representation of training species, if more than one training species exists in a cluster, cluster annotation is a concatention of names
### clusters with no training sequences in them are name c.[cluster number]
vicut_clusters <- read.table("454Reads.Lactobacillus.reduced.combined.filter.sing.hclust.minNodeCut.dir/minNodeCut.cltrs",header=T,as.is=T)
vicut_clusters$tax <- paste("c.",vicut_clusters$clstr,sep="")
vicut_clusters[which(!is.na(vicut_clusters$annot)),"tax"] <- paste(vicut_clusters[which(!is.na(vicut_clusters$annot)),"annot"],vicut_clusters[which(!is.na(vicut_clusters$annot)),"clstr"],sep=".")

annot_clusters <- tapply(vicut_clusters$annot[which(!is.na(vicut_clusters$annot))],vicut_clusters$clstr[which(!is.na(vicut_clusters$annot))],function(x) { y <- table(x); paste(names(y),collapse=".")})
annot_clusters <- data.frame(clusterID=as.numeric(names(annot_clusters)),Tax=paste(annot_clusters,names(annot_clusters),sep="."))

for(i in seq.int(1,nrow(annot_clusters))){
  vicut_clusters$tax[vicut_clusters$clstr == annot_clusters$clusterID[i]] <- annot_clusters$Tax[i]
}

cluster_mat$species <- vicut_clusters$tax[match(cluster_mat$X5,cluster_mat[na.exclude(match(vicut_clusters$readId,cluster_mat$X3)),"X5"])]

### Use FlashClust and WGCNA to produce clusters
library(flashClust)
library(WGCNA)
d5k <- read.table(dist.mat,skip=1,row.names=1)

hc <- flashClust(as.dist(d5k),method="average")

minModuleSize = 2;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = hc, distM = as.matrix(d5k),method="hybrid",
                            deepSplit = TRUE, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);


dynamicColors = labels2colors(dynamicMods)
#speciesColors = labels2colors(clusters$SPECIES)
#genusColors = labels2colors(clusters$GENUS)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#pdf("Dendrogram_fullsequence.pdf",width=36,height=8,pointsize=8)
plotDendroAndColors(hc, data.frame(TreeCut=dynamicColors), c("Tree Cut"),
                    dendroLabels = spp_taxon[match(hc$labels,spp_taxon[,1]),1], hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.dendroLabels=0.7,
                    main = "Clustering Lactobacillaceae 16S sequence")
#dev.off()


### Read in clusters computed by vicut, annotate by representation of training species, if more than one training species exists in a cluster, cluster annotation is a concatention of names
### clusters with no training sequences in them are name c.[cluster number]
sp <- spp_taxon[match(hc$labels,spp_taxon$ID),"Species"]
tax <- paste("c.",dynamicMods,sep="")
sp[which(is.na(sp))] <- tax[which(is.na(sp))]
clusterSpecies <- tapply(sp,tax,function(x) {tb <- table(x); if(length(tb) > 1) tb <- tb[-grep("c.",names(tb),fixed=T)];paste(names(tb),collapse=";")})
clusterSpecies <- clusterSpecies[match(tax,names(clusterSpecies))]

cluster_mat$species <- clusterSpecies[match(cluster_mat$X5,cluster_mat[na.exclude(match(hc$labels,cluster_mat$X3)),"X5"])]

## update the db, adding in cluster/species ID
update.species.rdp(con,cluster_mat)



