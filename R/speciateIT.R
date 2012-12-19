
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
project <- "Marmoset"
#project <- "JJ_Human_Vagina"
genus <- "Bifidobacterium"
nproc=16
mothur.template="/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta" 
output_dir="OutputFiles"
pipeline="dynamicTreeCut"

specieateMyReads <- function(project, genus, output_dir="OutputFiles", nproc = 8,pipeline=c("dynamicTreeCut","vicut"), mothur.template){

  pipeline <- match.arg(pipeline)
  ofile <- output.genus.reads(con,project=project,genus=genus, output_dir=output_dir)

## cluster using cdhit require 99.5% identity
  system(paste("cdhit-est -M 3000 -T 8 -d 0 -c 0.995 -n 9 -i ", ofile," -o ", gsub(".fasta",".reduced.0.995.fasta",ofile),sep=""))
  ## look at clustering results
  ofile <- gsub(".fasta",".reduced.0.995.fasta",ofile)
  cdhit_cluster <- readLines(paste(ofile,".clstr",sep=""))
  clust <- grep("^>Cluster",c(cdhit_cluster,">Cluster"))
  cdhit_cluster <- paste(cdhit_cluster,"C",rep(seq.int(1,length(diff(clust))),times=diff(clust)),sep=" ")
  cdhit_cluster <- cdhit_cluster[-clust]
  test_split <- strsplit(cdhit_cluster,split="\t|nt, >|\\.\\.\\. at +/|\\.\\.\\. |%C |C ")
  cluster_mat <- data.frame(matrix(unlist(test_split),ncol=5,byrow=T),stringsAsFactors=F)
  cluster_mat$X4 <- as.numeric(gsub("at [+]\\/|% ","",cluster_mat$X4))
  cluster_mat$errors <- as.numeric(cluster_mat$X2) - as.numeric(cluster_mat$X2)*(cluster_mat$X4/100)

  cluster_mat[is.na(cluster_mat$errors),"errors"] <- 0
  cluster_mat <- cluster_mat[order(as.numeric(cluster_mat$X5),-is.na(cluster_mat$X4)),]
 
########### PAWEL GAWER's SPECIATE IT PIPELINE -- MYVERSION
  ## add in Training species sequences
  ## maybe read them into R and then write out
  if (pipeline = "vicut"){  
    ifile <- ofile
    ofile <- gsub("0.995.fasta","combined.fa",ofile)
    speciateIThome <- "/mnt/home/msettles/opt/speciateit"
    speciesFile <- file.path(speciateIThome,"spp-data/Lactobacillus/rdp_vagi1200Len_Lactobacillus.V3V1.550bp_nr.fa")
    taxonFile <- file.path(speciateIThome,"spp-data/Lactobacillus/rdp_vagi1200Len_Lactobacillus.V3V1.550bp_nr.taxon")
    system(paste("cat ",speciesFile, ifile, ">", ofile,sep=" "))

    ## Use mothur to compute the Distance between sequences 
    #### Align the seuqences to silva
    #### Filter the alignment
    #### Compute distance matrix

    system(paste("mothur \"#align.seqs(candidate=",ofile,", template=", mothur.template ,", flip=T, processors=",nproc,"); ",
              "filter.seqs(fasta=",sub("fa","align",ofile),", processors=",nproc,"); ",
              "dist.seqs(fasta=",sub("fa","filter.fasta",ofile),", calc=onegap, output=square, processors=",nproc,");\"",sep=""))

    dist.mat <- sub("fa","filter.square.dist",ofile)
    outFile <- sub("fa","filter.sing.hclust.membStr",ofile)

    ############## Use FlashClust and Vicut to produce clusters
    library(flashClust)
    d5k <- read.table(dist.mat,skip=1,row.names=1)
    hc <- flashClust(as.dist(d5k),method="average")

    # using rowIds for leaves
    ### produces file ready for vicut
    hclust2merges2 <- function(hc,rowIds,filename)
    {
      sink(file=filename)
      internal_id = -1
      for (i in 1:nrow(hc$merge)){
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

    system(paste("vicut -t ", outFile, " -a ",taxonFile," -o ",sub("fa","filter.sing.hclust.minNodeCut.dir",ofile),sep=""))

    ### Read in clusters computed by vicut, annotate by representation of training species, if more than one training species exists in a cluster, cluster annotation is a concatention of names
    ### clusters with no training sequences in them are name c.[cluster number]
    vicut_clusters <- read.table(file.path(sub("fa","filter.sing.hclust.minNodeCut.dir",ofile),"minNodeCut.cltrs"),header=TRUE,as.is=TRUE)
    vicut_clusters$tax <- paste("c.",vicut_clusters$clstr,sep="")
    vicut_clusters[which(!is.na(vicut_clusters$annot)),"tax"] <- paste(vicut_clusters[which(!is.na(vicut_clusters$annot)),"annot"],vicut_clusters[which(!is.na(vicut_clusters$annot)),"clstr"],sep=".")

    annot_clusters <- tapply(vicut_clusters$annot[which(!is.na(vicut_clusters$annot))],vicut_clusters$clstr[which(!is.na(vicut_clusters$annot))],function(x) { y <- table(x); paste(names(y),collapse=".")})
    annot_clusters <- data.frame(clusterID=as.numeric(names(annot_clusters)),Tax=paste(annot_clusters,names(annot_clusters),sep="."))

    for(i in seq.int(1,nrow(annot_clusters))){
      vicut_clusters$tax[vicut_clusters$clstr == annot_clusters$clusterID[i]] <- annot_clusters$Tax[i]
    }

    cluster_mat$species <- vicut_clusters$tax[match(cluster_mat$X5,cluster_mat[na.exclude(match(vicut_clusters$readId,cluster_mat$X3)),"X5"])]
  }
############ Use FlashClust and WGCNA to produce clusters

  if (pipeline = "dynamicTreeCut"){
    system(paste("mothur \"#align.seqs(candidate=",ofile,", template=",mothur.template,", flip=T, processors=",nproc,"); filter.seqs(fasta=",gsub("fasta","align",ofile),", processors=",nproc,");\"",sep=""))
    system(paste("mv ",file.path(gsub("-","_",file.path(output_dir,project,"speciateIT")),"454Reads.filter"),file.path(gsub("-","_",file.path(output_dir,project,"speciateIT")),"454Reads.0.995.filter"),sep=" "))
    
    ## Prepare the reference
    filter <- readLines(file.path(gsub("-","_",file.path(output_dir,project,"speciateIT")),"454Reads.0.995.filter"))
    ref_seqs <- readBStringSet("Lactobacillaceae.patric.red.align") ## PreAligned
    keep <- which(unlist(strsplit(filter,"")[[1]]==1))
    #trimmed_ref <- DNAStringSet(gsub("[.]|-","",sapply(strsplit(as.character(ref_seqs),""),function(x) paste( x[keep] , collapse=""))))
    trimmed_ref <- BStringSet(sapply(strsplit(as.character(ref_seqs),""),function(x) paste( x[keep] , collapse="")))
    names(trimmed_ref) <-  paste("S",seq.int(1,length(trimmed_ref)),sep="")
    ifile <- gsub("fasta","filter.fasta",ofile)
    ofile <- gsub("filter.fasta","patric.filter.fasta",ifile)
    writeXStringSet(c(trimmed_ref,readBStringSet(ifile)),ofile)

    taxon_file <- read.table("Lactobacillaceae.patric.red.taxonomy",sep="\t",header=TRUE,as.is=TRUE)
    taxon_file <- data.frame(ID=paste("S",seq.int(1,length(trimmed_ref)),sep=""),Species=paste(taxon_file[,"GENUS"],taxon_file[,"SPECIES"],sep=" "))
    write.table(taxon_file,file=file.path(file.path(sub("-","_",output_dir,project,"speciateIT"),"Lactobacillaceae.patric.taxon")),sep="\t",row.names=FALSE,col.names=FALSE,quote=F)

    system(paste("mothur \"#dist.seqs(fasta=",ofile,", calc=onegap, output=square, processors=",nproc,");\"",sep=""))

    dist.mat <- sub("fasta","square.dist",ofile)

    library(flashClust)
    library(WGCNA)
    d5k <- read.table(dist.mat,skip=1,row.names=1)

    hc <- flashClust(as.dist(d5k),method="average")

    minModuleSize = 1;
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = hc, distM = as.matrix(d5k),method="hybrid",
                            deepSplit = TRUE, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

    dynamicColors = labels2colors(dynamicMods)

#    table(dynamicColors)
    # Plot the dendrogram and colors underneath
    pdf(file.path(sub("-","_",file.path(output_dir,project,"speciateIT")),"Dendrogram_dynamicTreeCut.pdf"),width=36,height=8,pointsize=8)
    plotDendroAndColors(hc, data.frame(TreeCut=dynamicColors), c("Dynamic Tree Cut"),
                    dendroLabels = taxon_file[match(hc$labels,taxon_file[,1]),1], hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.dendroLabels=0.7,
                    main = paste("Clustering",genus,"16S sequence",sep=""))
    dev.off()

    sp <- taxon_file[match(hc$labels,taxon_file$ID),"Species"]
    tax <- paste("c.",dynamicMods,sep="")
    sp[which(is.na(sp))] <- tax[which(is.na(sp))]
    clusterSpecies <- tapply(sp,tax,
                         function(x) {
                           tb <- table(x); 
                           if(length(tb) > 1) {
                             rm <- grep("c.",names(tb),fixed=T);
                             if (length(rm)) tb <- tb[-rm];
                           }
                           paste(names(tb),collapse=";")}
                         )
    clusterSpecies <- clusterSpecies[match(tax,names(clusterSpecies))]

    cluster_mat$species <- clusterSpecies[match(cluster_mat$X5,cluster_mat[match(hc$labels,cluster_mat$X3),"X5"])]
  }
  ## update the db, adding in cluster/species ID
  update.species.rdp(con,cluster_mat)
}
########################################################################################################
########################################################################################################
### PLOTTING DIFFERENCES BETWEEN VICUT AND DYNAMICTREECUT

sql <- paste("SELECT read_data.*, rdp_report.*
                             FROM rdp_report, read_data
                             WHERE read_data.Acc IN ($X3)
             AND rdp_report.lucyUnique = read_data.lucyUnique",sep="")

rdp <- dbGetPreparedQuery(con,sql,bind.data=cluster_mat)
rdp2 <- rdp[match(hc$labels,rdp$Acc),]
cluster_mat2 <- cluster_mat[match(hc$labels,cluster_mat$X3),]

vicutColors <- labels2colors(cluster_mat2$species_vicut_single)
dynamicColors <- labels2colors(cluster_mat2$species)

rwgPalette <- redWhiteGreen(100, gamma = 1)

#This adds a column of color values
# based on the y values
LucyLength <- rwgPalette[as.numeric(cut(rdp2$LucyLength,breaks = 100))]
LucyNs <- labels2colors(rdp2$LucyNs)
LucymHomoPrun <- labels2colors(rdp2$LucymHomoPrun)
genus_bootstrap <- rwgPalette[as.numeric(cut(rdp2$genus_bootstrap,breaks=100))]

pdf("Dendrogram_Comparisons.pdf",width=36,height=8,pointsize=8)

plotDendroAndColors(hc, data.frame(TreeCut=dynamicColors,ViCut=vicutColors,LucyLength,LucyNs,LucymHomoPrun,genus_bootstrap), c("Dynamic Tree Cut","ViCut","Read Length","Ns","Homopolymer run","genus_bootstrap"),
                    dendroLabels = taxon_file[match(hc$labels,taxon_file[,1]),1], hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.dendroLabels=0.7,
                    main = paste("Clustering",genus,"16S sequence",sep=""))
dev.off

d5l <- d5k[-(grep("^S",rownames(d5k))),-grep("^S",rownames(d5k))]
mds <- cmdscale(as.dist(d5l), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

pdf("mds-dynamicTree.pdf",width=9,height=9,pointsize=8)
plot(mds[,1], mds[,2], type = "p",col=dynamicColors, pch=20, xlab = "", ylab = "", asp = 1, axes = FALSE,
     main = "cmdscale(16S lactobacillaea")
dev.off()

pdf("mds-vicut.pdf",width=9,height=9,pointsize=8)
plot(mds[,1], mds[,2], type = "p",col=vicutColors, pch=20, xlab = "", ylab = "", asp = 1, axes = FALSE,
     main = "cmdscale(16S lactobacillaea")
dev.off()

### See if we can further deliniate sequences and check

alignment <- readBStringSet("OutputFiles/Witkin_VVS/speciateIT/454Reads.Lactobacillus.reduced.0.995.patric.filter.fasta")

##########################################################################################
### TEST C-VALUES
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


