## load libraries
library(Biobase)
library(RColorBrewer)
library(vegan)
library(WGCNA)
library(vegetarian)
library(gplots)

name <- "Witkin-VVS"
source("heatmap.microbe.R")

heatcol <- rainbow(256,start=0.2,end=1)
## SET up data
## Read in, remove unnecessary columns, etc.
abund2 <- read.table(paste(name,"abundance","txt",sep="."),sep="\t",header=T,row.names=1)

da <- dim(abund2)
abund2 <- abund2[-da[1],] ### remove column totals
abund2 <- abund2[,-da[2]] #### romove row totals
level <- abund2$Level   ### set aside Levels
abund2$Level <- NULL ### Remove

metaData <- read.table("Jacobson_MetaData.ordered.txt",sep=" ",header=TRUE)

readData <- read.table("JJ_Human_Vagina.readcounts.txt",sep="\t",header=TRUE)

metaData$fail[match(levels(factor(readData$Sample_ID)),metaData$Sample_ID)] <- tapply(readData$Fail,readData$Sample_ID,sum)
metaData$pass[match(levels(factor(readData$Sample_ID)),metaData$Sample_ID)] <- tapply(readData$Pass,readData$Sample_ID,sum)

mmatch <- match(paste("X",metaData$Sample_ID,sep=""),colnames(abund))
abund <- abund[,match] ### ensure correct order

rem <- which(rowSums(abund)==0)
if (length(rem) > 0) {abund <- abund[-rem,]; level <- level[-rem];}

### filter Samples
metaData$Read_Count <- colSums(abund)
remSample <- colSums(abund) >= 1000 & !is.na(colSums(abund))
abund <- abund[,remSample]
metaData <- metaData[remSample,]

metaData$read_count <- colSums(abund)
colnames(abund) <- metaData$Sample_ID
#abund <- abund[order(rowSums(abund),decreasing=T),]
freqAbund <- sweep(abund2,2,colSums(abund2,na.rm=TRUE),"/")
#########################################################################################################################
#### WRITE OUT NEW DATA
### compute frequencies and write out

## filter out genus for table
keepRows <- (apply(freqAbund >= 0.01,1,sum,na.rm=TRUE) > 1) | (apply(freqAbund >= 0.05, 1, sum,na.rm=TRUE) > 0) ## remove taxa with < 0.05 and only found in 1
keepRows["Bacteria"] <- FALSE

Other <- colSums(abund2[!keepRows,])
abund <- abund2[keepRows,]
level <- level[keepRows]

abund <- rbind(abund,Other=Other)
level <- c(level,"Bacteria")
freqAbund <- sweep(abund,2,colSums(abund),"/")

write.table(data.frame(Names=rownames(freqAbund),Level=level,freqAbund),paste(name,"analysis.relative.abundance.REDUCED",sep="."),sep="\t",row.names=FALSE,col.names=TRUE)
write.table(data.frame(Names=rownames(abund),Level=level,abund),paste(name,"analysis.abundance",sep="."),sep="\t",row.names=FALSE,col.names=TRUE)
write.table(metaData,file=paste(name,"analysis.metaData.txt",sep="."),sep="\t",row.names=TRUE,col.names=TRUE)

md_col <- c("Sample_ID","Subject","Visit","Source","Replicate")
write.table(t(metaData[,md_col]),file=paste(name,"analysis.relative.full.abundance",sep="."),sep="\t",row.names=TRUE,col.names=FALSE)
write.table(data.frame(Names=rownames(freqAbund),Level=level,round(freqAbund,digits=2)),paste(name,"analysis.relative.full.abundance",sep="."),sep="\t",row.names=FALSE,col.names=TRUE,append=TRUE)


write.table(t(metaData[,md_col]),file=paste(name,"analysis.full.abundance",sep="."),sep="\t",row.names=TRUE,col.names=FALSE)
write.table(data.frame(Names=rownames(abund),Level=level,abund),paste(name,"analysis.full.abundance",sep="."),sep="\t",row.names=FALSE,col.names=TRUE,append=TRUE)

#########################################################################################################################


pp <- ((rowSums(abund,na.rm=T))/sum(abund,na.rm=T))
write.table(data.frame(Names=rownames(freqAbund[pp>=0.01,]),Level=level[pp>=0.01],freqAbund[pp>=0.01,]),paste(name,"analysis.relative.1percent.abundance",sep="."),sep="\t",row.names=FALSE,col.names=TRUE)
write.table(data.frame(Names=rownames(freqAbund[pp>=0.05,]),Level=level[pp>=0.05],freqAbund[pp>=0.05,]),paste(name,"analysis.relative.5percent.abundance",sep="."),sep="\t",row.names=FALSE,col.names=TRUE)




Full_Meta_Data <- metaData
library(plotrix)
library(gplots)

my_colors <- unique(sub("[0-9]+$","",colors()))
load("phColorTbl_ct.2k.RData")
my_colors <- setdiff(my_colors,unique(sub("[0-9]+$","",phColorTbl)))
names(phColorTbl) <- sub("^L_","Lactobacillus ",names(phColorTbl))
names(phColorTbl) <- sub("_[0-9]+$","",names(phColorTbl))

taxa_colors <- phColorTbl[match(rownames(abund),names(phColorTbl))]
taxa_colors[is.na(taxa_colors)] <- my_colors[sample(1:length(my_colors),sum(is.na(taxa_colors)))]
names(taxa_colors) <- rownames(abund)

pdf(file=paste("sampleplots.new.pdf"),width=6,height=9,pointsize=8)
for (i in sort(unique(Full_Meta_Data$Subject))){
  #i = "2"
  sample <- abund[,which(Full_Meta_Data$Subject == i)]
  freq <- sweep(sample,2,colSums(sample),"/")
  meta <- Full_Meta_Data[which(Full_Meta_Data$Subject == i),]
  ord <- order(meta$Source,meta$Visit)
  freq <- freq[,ord]
  meta <- meta[ord,]

  layout(matrix(c(1,1,1,4,1,1,1,4,2,2,2,4,2,2,2,4,3,3,3,4,3,3,3,4,0,0,0,0),nrow=7,byrow=T))

  par(mar=c(2, 4, 3, 2) + 0.1)
  Top20 <- order(rowMeans(as.matrix(sample[,meta$Source=="Vagina"])),decreasing=T)[1:10]
  fullTaxa <- Top20
  stackpoly(t(freq[Top20,meta$Source=="Vagina"]),stack=T,col=taxa_colors[Top20],xaxlab=meta$Visit[meta$Source=="Vagina"], ylim=c(0,1),main="Vagina",ylab="Proportion" )
  axis(meta$Replicate[meta$Source=="Vagina"],at=1:length(meta$Replicate[meta$Source=="Vagina"]),side=1,line=1,tick=F)

  Top20 <- order(rowMeans(as.matrix(sample[,meta$Source=="Cervix"])),decreasing=T)[1:10]
  fullTaxa <- union(fullTaxa,Top20)
  stackpoly(t(freq[Top20,meta$Source=="Cervix"]),bg="grey",stack=T,col=taxa_colors[Top20],xaxlab=meta$Visit[meta$Source=="Cervix"],ylim=c(0,1),main="Cervix",ylab="Proportion")
  axis(meta$Replicate[meta$Source=="Cervix"],at=1:length(meta$Replicate[meta$Source=="Cervix"]),side=1,line=1,tick=F)

  Top20 <- order(rowMeans(as.matrix(sample[,meta$Source=="Uterine"])),decreasing=T)[1:10]
  fullTaxa <- union(fullTaxa,Top20)
  stackpoly(t(freq[Top20,meta$Source=="Uterine"]),stack=T,col=taxa_colors[Top20],xaxlab=meta$Visit[meta$Source=="Uterine"],ylim=c(0,1),main="Uterine",xlab="Treatment",ylab="Proportion")
  axis(meta$Replicate[meta$Source=="Uterine"],at=1:length(meta$Replicate[meta$Source=="Uterine"]),side=1,line=1,tick=F)
  mtext("Visit",side=1,line=3,cex=0.8)
  mtext(paste("Subject",i),side=1,line=5)
  par(mar = c(5,0,4,0))

  ord <- order(rownames(freq)[fullTaxa])
  plot(seq(1,10), seq(1,10), type = "n", axes = F,xlab="",ylab="")
  legend("bottomleft", rownames(freq)[fullTaxa][ord],bty="n",bg="#ffffff55",fill=taxa_colors[fullTaxa][ord],inset=0)

#  QFM <- meta[1,c(2,,8,9,10,11)]
#  textplot("Subject 2",show.rownames=FALSE,show.colnames = FALSE,valign="top")

}
dev.off()


freqa <- sweep(abund,2,colSums(abund),"/")
pdf("TaxaRepresentation.pdf",width=7,height=10,pointsize=8)
par(mar=c(3,10,4,2))
for (i in unique(Full_Meta_Data$Source)){
   redmat <- freqa[,Full_Meta_Data$Source == i]
   redmeta <- Full_Meta_Data[Full_Meta_Data$Source == i,]
   maintext <- paste("Taxa represention across sample (",redmeta[1,"Source"],")")
   boxplot(t(redmat[order(rowSums(redmat),decreasing=F),]),horizontal=T,las=1,cex=0.3,cex.axis=0.75,
main=maintext)
}

dev.off()

### done preparing data
##################################################################################################


pdf(file="cluster-analysis.pdf",width=18,height=5.5,pointsize=6)
par(oma=c(2,2,2,2))
for (i in unique(Full_Meta_Data$Source)){
  redmat <- freqa[,Full_Meta_Data$Source == i]
  abdist <- vegdist(t(redmat),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
  redmeta <- Full_Meta_Data[Full_Meta_Data$Source == i,]
  maintext <- paste("Cluster Dentrogram (",redmeta[1,"Source"],")")
  hcl <- hclust(abdist,method="average")
  plot(hcl,labels=redmeta[,"Sample_ID"],family="sans",cex=0.75,font=2,main=maintext,sub="Bray's clustering, average linkage" ,ylab="Distance",xlab="")
}
dev.off()

pdf(file="heatmap-analysis.pdf",width=18,height=10,pointsize=10)
par(oma=c(2,2,2,2))
for (i in unique(Full_Meta_Data$Source)){
  redmat <- freqa[,Full_Meta_Data$Source == i]
  abdist <- vegdist(t(redmat),method="bray",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
  redmeta <- Full_Meta_Data[Full_Meta_Data$Source == i,]
  maintext <- paste("Cluster Dentrogram (",redmeta[1,"Source"],")")
  hcl <- hclust(abdist,method="average")
  
  rsum <- apply(redmat,1,max)
  pick <- order(rsum,decreasing=TRUE)[1:25]
  
  par(oma=c(2,2,4,10))
  heatmap.microbe(as.matrix(redmat)[pick,hcl$order],
                  col=heatcol,
                  distfun=vegdist,
                  labCol=subjectstring,
                  #                ColSideColors=fcol,
                  #                ColSideLabels="Cervix/Vagina",
                  #                colsep=which(!duplicated(mdata$Subject)[-1]),
                  cexRow=1.3,
                  cexCol=0.8,
                  mar=c(7,5),
                  dendrogram="none",
                  Rowv=F,
                  Colv=F,
                  sepcolor="white",
                  keysize=0.5,
                  trace="none",
                  key=T,
                  density.info="none",family="sans")
  
  mtext(paste("Heatmap of top 25 taxa present (",redmeta[1,"Source"],")"),side=3,line=1)  
}
dev.off()

