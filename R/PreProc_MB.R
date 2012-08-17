#! /usr/bin/env Rscript --vanilla --default-packages=utils,rSFFreader --verbose
###############################################################################
###
### Primary functions/main for Amplicon Processing
###
###############################################################################
options(width=120)
options("stringsAsFactors" = FALSE)
library(rSFFreader)
library(getopt)

source("functions.R")

#sfffiles <- commandArgs(TRUE)

sfffiles <- "../Amplicon_SFFfiles_AmpProc/GKF1FGD01.sff"
screenfile <- "screen_V3-V1Jun11.fa"
tagkey <- "^FP_"
tag_nucs <- 8
primerkey <- "^RP_"

maxNs <- 2
maxforwardprimererrors <- 2
maxhammingdisttag <- 1
minlength <- 100
maxlength <- 600
maxhomopol <- 7

rdpPath <- "/mnt/home/msettles/opt/rdp_classifier_2.2/rdp_classifier-2.2.jar"
nproc=4
## for speciation
reverseSeq <- TRUE
mothur.template="/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta" 


#### This script preprocesses Roche 454 Reads
## START Rdev

basefilename <- sub(".sff","",basename(sfffiles))
checksystem()

outfile <- paste(basefilename,".out",sep="")


fq <- readsff(sfffiles)
cat(paste("Total Number of Reads in SFF files:",length(fq),"\n"),file=outfile,append=F)    
clipMode(fq) <- "Raw"

#####################
### Base Roche stats and clip
print("Identifying Roche Clip points")
ReadData <- data.frame(Acc=as.character(id(fq)),
                       Run=substring(ReadData$Acc,1,9),
                       RawLength=width(fq))

ReadData$RocheLC <- start(qualityClip(fq))
ReadData$RocheRC <- end(qualityClip(fq))

ReadData$RocheLength <- width(sread(fq,clipmode="Full"))

cat(paste("Mean length before Roche Right Clip",signif(mean(ReadData$RawLength),3),"\n"),file=outfile,append=T)    
cat(paste("Median length before Roche Right Clip",signif(median(ReadData$RawLength),3),"\n"),file=outfile,append=T)    
cat(paste("Mean length after Roche Right Clip",signif(mean(ReadData$RocheLength),3),"\n"),file=outfile,append=T)    
cat(paste("Median length after Roche Right Clip",signif(median(ReadData$RocheLength),3),"\n"),file=outfile,append=T)    

### End Base Roche
#####################
###
#####################
### Run cross_match looking for primers and tags

writeXStringSet(sread(fq), "TMP.raw.fasta",append=FALSE,format="fasta")
writePhredQual(fq,"TMP.raw.fasta.qual",mode="w")

#screenfile  <- "screen_V3-V1Jun.fa"
### USE minscore >15 need to test <= 15
system(paste("cross_match TMP.raw.fasta ../", screenfile, " -minmatch 12 -minscore 11 -tags > TMP.cmout",sep=""))
cm_out <- parse_cm("TMP.cmout")

cm_outLess14 <- cm_out[which(cm_out$score <= 14),]
cm_out <- cm_out[-which(cm_out$score <= 14),]
### End CrossMatch
############################
###
############################
### Adapter Clip Points
#ReadData$AdapterLC <- ReadData$RocheLC
#ReadData$AdapterRC <- ReadData$RocheRC
ReadData$AdapterLC <- 5
ReadData$AdapterRC <- ReadData$RawLength

print("Parsing cross_match output")
### SET new clip points, recording primer information

#### New Left Clip point starts after the inner most Forward Primer

cm_out$read_end[(cm_out$FC == "F")] = cm_out$read_end[(cm_out$FC == "F")] + cm_out$adapt_remain[(cm_out$FC == "F")] + 1

tmp <-tapply(cm_out$read_end[(cm_out$FC == "F")],cm_out$read_id[(cm_out$FC == "F")],max)
ReadData$AdapterLC[match(names(tmp),ReadData$Acc)] <-
  pmax(ReadData$AdapterLC[match(names(tmp),ReadData$Acc)],tmp)

ReadData$AdapterLC <- pmin(ReadData$RawLength,ReadData$AdapterLC)

#### New Right Clip points occurs before the inner most Reverse Primers
cm_out$read_start[(cm_out$FC == "C")] = cm_out$read_start[(cm_out$FC == "C")] - cm_out$adapt_remain[(cm_out$FC == "C")] - 1
cm_out$read_start[cm_out$read_start <=0] <- 5
tmp <-tapply(cm_out$read_start[(cm_out$FC == "C")],cm_out$read_id[(cm_out$FC == "C")],min)
ReadData$AdapterRC[match(names(tmp),ReadData$Acc)] <-
  pmin(ReadData$AdapterRC[match(names(tmp),ReadData$Acc)],tmp)

ReadData$AdapterRC <- pmin(ReadData$RawLength,ReadData$AdapterRC)
ReadData$AdapterLC <- pmin(ReadData$AdapterLC,ReadData$AdapterRC)
## idendify reads with duplicate matches to Forward or Reverse Primers
readsWdups <- unique(c( cm_out$read_id[which(cm_out$FC == "F")][duplicated(cm_out$read_id[which(cm_out$FC == "F")])], ## forward
                        cm_out$read_id[which(cm_out$FC == "C")][duplicated(cm_out$read_id[which(cm_out$FC == "C")])])) ## reverse

cat(paste("Identified ::", length(readsWdups),":: reads with multiple adapters or primers, these are likely primer dimer sequences and are removed\n"),file=outfile,append=T)


## remove reads that have multiple matches to forward or reverse probes
if (length(readsWdups) > 0){
  cm_out <- cm_out[-which(cm_out$read_id %in% readsWdups),]
  ReadData$AdapterRC[match(readsWdups,ReadData$Acc)] <- ReadData$AdapterLC[match(readsWdups,ReadData$Acc)]
}
## remove all reads where forward primer does not start at base 1
## requiring that the tagged primer start at 1
cat(paste("removing Forward matches that do not start at base 1 :",length(which(cm_out$read_start != 1 & cm_out$FC == "F")), " reads removed\n"),file=outfile,append=T)
### fix, remove reads##
cm_out <- cm_out[-which(cm_out$read_start != 1 & cm_out$FC == "F"),]

ReadData$AdapterLength <- ReadData$AdapterRC - ReadData$AdapterLC +1

clipMode(fq) <- "Raw"
##########################
## Calculate Primer errors INCLUDES TAG!!
print("Computing, tags and tagged primer errors")
primer_offset <-  as.numeric(names(sort(table(cm_out$adapt_start[which(cm_out$FC == "F")]),decreasing=T))[1])
ReadData$Primer_Code <- NA
ReadData$Primer_Code[match(cm_out$read_id[grep(tagkey,cm_out$adapt)],ReadData$Acc)] <- cm_out$adapt[grep(tagkey,cm_out$adapt)]
### ADD IN DUMMY BARCODE
ReadData$Barcode <- ReadData$Primer_Code
ReadData$FPErr <- 100
ReadData$FPErr[match(cm_out$read_id[grep(tagkey,cm_out$adapt)],ReadData$Acc)] <- 
  cm_out[grep(tagkey,cm_out$adapt),"err"] - primer_offset +1

cat("Forward Primer Errors:\n",file=outfile,append=T)
sink(file=outfile,append=TRUE)
stem(ReadData$FPErr,scale=0.5,width=10)
sink()
cat("\n",file=outfile,append=T)

print("Computing, primer tag hamming distance from assigned tag")
screen <- read.DNAStringSet(file.path("..",screenfile))
if(screenfile == "screen_V3-V1Jun11.fa"){
  tag_nucs <- as.numeric(sapply(strsplit(names(screen),split=","),function(x) x[2]))
  tags <- subseq(screen,primer_offset+4,width=tag_nucs)[match(ReadData$Primer_Code[!is.na(ReadData$Primer_Code)],names(screen))]
  rtags <- subseq(sread(fq)[!is.na(ReadData$Primer_Code)],5,width=tag_nucs[match(ReadData$Primer_Code[!is.na(ReadData$Primer_Code)],names(screen))])
  tagpairs <- cbind(as.character(tags),as.character(rtags))
  tagdiff <- numeric(nrow(tagpairs))  ## Sets all to 0
  tagdiff[tagpairs[,1] != tagpairs[,2]] <- apply(tagpairs[tagpairs[,1] != tagpairs[,2],],1,stringDist, method="hamming")
  ReadData$Code_Dist[!is.na(ReadData$Primer_Code)] <- tagdiff
}else{
  tags <- subseq(screen,primer_offset+4,width=tag_nucs)[match(ReadData$Primer_Code[!is.na(ReadData$Primer_Code)],names(screen))]
  rtags <- subseq(sread(fq)[!is.na(ReadData$Primer_Code)],5,width=tag_nucs)
  tagpairs <- cbind(as.character(tags),as.character(rtags))
  tagdiff <- numeric(nrow(tagpairs))  ## Sets all to 0
  tagdiff[tagpairs[,1] != tagpairs[,2]] <- apply(tagpairs[tagpairs[,1] != tagpairs[,2],],1,stringDist, method="hamming")
  ReadData$Code_Dist[!is.na(ReadData$Primer_Code)] <- tagdiff
}

cat("Tag Hamming Distance:\n",file=outfile,append=T)
sink(file=outfile,append=TRUE)
stem(ReadData$Code_Dist,scale=0.5,width=10)
sink()
cat("\n",file=outfile,append=T)

print("Determining end adapter presence and error")
ReadData$Primer_Reverse <- NA
ReadData$Primer_Reverse[match(cm_out$read_id[grep(primerkey,cm_out$adapt)],ReadData$Acc)] <- cm_out$adapt[grep(primerkey,cm_out$adapt)]
ReadData$RPErr <- 100
ReadData$RPErr[match(cm_out$read_id[grep(primerkey,cm_out$adapt)],ReadData$Acc)] <- 
  cm_out[grep(primerkey,cm_out$adapt),"err"]

cat("Reverse Primer Distance\n",file=outfile,append=T)
sink(file=outfile,append=TRUE)
stem(ReadData$RPErr,scale=0.5,width=10)
sink()
cat("\n",file=outfile,append=T)

cat(paste("Adapter Based Clipping\n"),file=outfile,append=T)
cat(paste("Total Forward Primer Errors <= :",maxforwardprimererrors, " :: Reads meeting that criteria ::",table(ReadData$FPErr <= maxforwardprimererrors)[2], "\n"),file=outfile,append=T)
cat(paste("Tag Max Hamming Distance from Target <= :",maxhammingdisttag, " :: Reads meeting that criteria ::",table(ReadData$Code_Dist <= maxhammingdisttag)[2], "\n"),file=outfile,append=T)
#cat(paste("Maximum Number of Abiguous Characters (Ns) <= :",maxNs, " :: Reads meeting that criteria ::",table(ReadData$AdapterNs <= maxNs)[2], "\n"),file=outfile,append=T)
#cat(paste("Maximum Length of a Homopolymer Run <= :",maxhomopol, " :: Reads meeting that criteria ::",table(ReadData$mHomoPrun <= maxhomopol)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Greater than Min Length > :",minlength, " :: Reads meeting that criteria ::",table(ReadData$AdapterLength > minlength)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Less than Max Length :",maxlength, " :: Reads meeting that criteria ::",table(ReadData$AdapterLength < maxlength)[2], "\n"),file=outfile,append=T)
cat(paste("Primer with Tab found :",table(!is.na(ReadData$Primer_Code))[2], "\n"),file=outfile,append=T)

#####################################################
ReadData$keepAdapter <- FALSE
ReadData$keepAdapter <- #keep[expand] &
  ReadData$FPErr <= maxforwardprimererrors & 
  ReadData$Code_Dist <= maxhammingdisttag &
  #                        ReadData$AdapterNs <= maxNs &
  #                        ReadData$mHomoPrun <= maxhomopol & 
  ReadData$AdapterLength > minlength &
  ReadData$AdapterLength < maxlength &
  !is.na(ReadData$Primer_Code)

ReadData$keepAdapter[is.na(ReadData$keepAdapter)] <- FALSE
cat(paste("RESULTS\nAdapter based filter: ",sum(ReadData$keepAdapter), "Reads kept ::", sum(!ReadData$keepAdapter), "Reads removed\n"),file=outfile,append=T)


### Write out data
primerRange <- IRanges(start=ReadData$AdapterLC,ReadData$AdapterRC)
adapterClip(fq) <- primerRange ## need to add replacement function 
fq@adapterClip <- primerRange
clipMode(fq)  <- "Full"

writeXStringSet(sread(fq), "TMP.adapterClip.fasta",append=FALSE,format="fasta")
writePhredQual(fq,"TMP.adapterClip.fasta.qual",mode="w")

###########################
## check for poor regions with lucy
print("Performing lucy based quality trimming, error (27)")
system("lucy -xtra 4 -debug TMP.lucy_clip.txt -error 0.002 0.002 -output TMP.lucy.fasta TMP.lucy.fasta.qual  TMP.adapterClip.fasta  TMP.adapterClip.fasta.qual")

lucy <- read.table("TMP.lucy_clip.txt",as.is=T)

ReadData$lucyLC <- ReadData$AdapterLC
ReadData$lucyRC <- ReadData$AdapterRC

zeros <- which(lucy$V9 == 0)
ReadData$lucyRC[match(lucy[zeros,1] ,ReadData$Acc)] <- ReadData$lucyLC[match(lucy[zeros,1] ,ReadData$Acc)]
ReadData$lucyRC[match(lucy[-zeros,1],ReadData$Acc)] <- ReadData$AdapterLC[match(lucy[-zeros,1],ReadData$Acc)] + lucy[-zeros,10] -1
ReadData$lucyLC[match(lucy[-zeros,1],ReadData$Acc)] <- ReadData$AdapterLC[match(lucy[-zeros,1],ReadData$Acc)] + lucy[-zeros,9] -1
ReadData$lucyLength <- ReadData$lucyRC - ReadData$lucyLC +1

print("computer read statistics")
lucyRange <- IRanges(start=ReadData$lucyLC,ReadData$lucyRC)
fq@adapterClip <- lucyRange
clipMode(fq)  <- "Full"

fa.lucy <- sread(fq)
## get unique sequences
ReadData$LucyUnique <- NA
fa.lucy <- sort(fa.lucy)
lucy.unique <- length(unique(fa.lucy))
UniqueID <- rep(paste(basefilename,seq.int(lucy.unique),sep="."),times=diff(c(which(!duplicated(fa.lucy)),length(fa.lucy)+1)))
ReadData$LucyUnique[match(names(fa.lucy),ReadData$Acc)] <- UniqueID

cat(paste("total number of unique reads:", lucy.unique,"\n"),file=outfile,append=T)

expand <- match(ReadData$LucyUnique,unique(UniqueID))

lucyAlphFreq <- alphabetFrequency(fa.lucy[!duplicated(UniqueID)])

ReadData$lucyNs <- lucyAlphFreq[expand,"N"]
homoRuns <-  sapply(fa.lucy[!duplicated(UniqueID)],function(x) sapply(c("A","C","T","G","N"),function(y) longestConsecutive(as.character(x),letter=y)))
ReadData$lucymHomoPrun <- apply(t(homoRuns),1,max)[expand]

cat(paste("\nLucy Additional Clipping\n"),file=outfile,append=T)
cat(paste("Maximum Number of Abiguous Characters (Ns) <= :",maxNs, " :: Reads meeting that criteria ::",table(ReadData$lucyNs <= maxNs)[2], "\n"),file=outfile,append=T)
cat(paste("Maximum Length of a Homopolymer Run <= :",maxhomopol, " :: Reads meeting that criteria ::",table(ReadData$lucymHomoPrun <= maxhomopol)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Greater than Min Length > :",minlength, " :: Reads meeting that criteria ::",table(ReadData$lucyLength > minlength)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Less than Max Length :",maxlength, " :: Reads meeting that criteria ::",table(ReadData$lucyLength < maxlength)[2], "\n"),file=outfile,append=T)


########
## Run the RDP classifier on all the data
######### Lucy Clipped

fq_uni <- fq[(!duplicated(ReadData$LucyUnique)& ReadData$lucyLength > 50)]
#fq_sread <- sread(fq_uni)
#names(fq_sread) <- ReadData$LucyUnique[(!duplicated(ReadData$LucyUnique)& ReadData$lucyLength > 50)]
#fq_uni <- SffReadsQ(fq_sread,quality(fq_uni), qualityClip(fq_uni), adapterClip(fq_uni), clipMode = "Full")
writeXStringSet(sread(fq_uni), "TMP.rdp.fasta",append=FALSE,format="fasta")
writePhredQual(fq_uni,"TMP.rdp.fasta.qual",mode="w")

######################
system(paste("java -Xmx1g -jar ",rdpPath," -q TMP.rdp.fasta -o TMP.lucy.rdpV6.fix -f fixrank",sep=""))
rdp.lucy <- read.table("TMP.lucy.rdpV6.fix",sep="\t")
expand <- match(ReadData$LucyUnique,ReadData$LucyUnique[match(rdp.lucy[,1],ReadData$Acc)])
ReadData$LucyRDPgenus <- rdp.lucy[expand,"V21"]
ReadData$LucyRDPboot <- rdp.lucy[expand,"V23"]

system(paste("mothur \"#align.seqs(candidate=TMP.rdp.fasta, template=", mothur.template ,", flip=T, processors=",nproc,")\"",sep=""))
ReadData$LucyFlip = FALSE
flip <- read.table("TMP.rdp.flip.accnos",sep="\t")

if(ncol(flip) > 1){
  expand <- match(ReadData$LucyUnique,
                  ReadData$LucyUnique[match(flip[grep("reverse complement produced a better alignment",flip[,2]),1],ReadData$Acc)])
  ReadData$LucyFlip <- !is.na(expand)
}

align.report <- read.table("TMP.rdp.align.report",sep="\t",header=T)
expand <- match(ReadData$LucyUnique,ReadData$LucyUnique[match(align.report[,1],ReadData$Acc)])
align.report$QueryFull <- ReadData$lucyLength[match(align.report[,1],ReadData$Acc)]
align.report$flip <- ReadData$LucyFlip[match(align.report[,1],ReadData$Acc)]

align.report$adpLC <- ReadData$AdapterLC[match(align.report[,1],ReadData$Acc)]


keep <- align.report$QueryStart < 5 &
  (align.report$QueryFull-align.report$QueryLength) < 5 &
  align.report$flip == reverseSeq
#            align.report$QueryLength > 200 &    ### not sure if we need this

#TEtm <-  mean(align.report$TemplateEnd[keep],trim=0.1)
#TEsd <-  sd.trim(align.report$TemplateEnd[keep],trim=0.1)

#cat(paste("Lucy Clipped Alignment Statistics\n"),file=outfile,append=T)
#cat(paste("Template End mean:",signif(TEtm,3),":Template End Sd:",signif(TEsd,3),"\n"),file=outfile,append=T)    

#keep <- keep &
#            align.report$TemplateEnd > (TEtm - 4*TEsd) & 
#            align.report$TemplateEnd < (TEtm + 4*TEsd)

TEtm <- 510 ### set static align end position
keep <- keep &
  align.report$TemplateEnd > (TEtm - 75) & 
  align.report$TemplateEnd < (TEtm + 75)


ReadData$LucyTE <- align.report$TemplateEnd[expand]
ReadData$LucyQM <- align.report$QueryLength[expand]


ReadData$keepLucy <- FALSE
ReadData$keepLucy[match(lucy[,1],ReadData$Acc)] <- TRUE
ReadData$keepLucy <-    keep[expand] &
  ReadData$FPErr <= maxforwardprimererrors &
  ReadData$Code_Dist <= maxhammingdisttag &
  ReadData$lucyNs <= maxNs &
  ReadData$lucymHomoPrun <= maxhomopol &
  ReadData$lucyLength > minlength &
  ReadData$lucyLength < maxlength &
  !is.na(ReadData$Primer_Code)

ReadData$keepLucy[is.na(ReadData$keepLucy)] <- FALSE

cat(paste("RESULTS\nAdditional Lucy clipped filter: ",sum(ReadData$keepLucy), "Reads kept ::", sum(!ReadData$keepLucy), "Reads removed\n"),file=outfile,append=T)

###########################
cat("Overlap between Adapter Clipped and Lucy Clipped Data:\n",file=outfile,append=T)
sink(file=outfile,append=TRUE)
table(Adapter=ReadData$keepAdapter,Lucy=ReadData$keepLucy)
sink()
cat("\n",file=outfile,append=T)
cat(paste("Mean number of bases removed by Lucy clipping:",signif(mean(ReadData$AdapterLength-ReadData$lucyLength),3),"\n"),file=outfile,append=T)
cat(paste("Median number of bases removed by Lucy clipping:",signif(median(ReadData$AdapterLength-ReadData$lucyLength),3),"\n"),file=outfile,append=T)

cat(paste("Number of reads Adapter Clipped and Unique:",table(!duplicated(ReadData$AdapterUnique) & ReadData$keepAdapter)[2],"\n"),file=outfile,append=T)
cat(paste("Number of reads Lucy    Clipped and Unique:",table(!duplicated(ReadData$LucyUnique) & ReadData$keepLucy)[2],"\n"),file=outfile,append=T)
##########################

ReadData$Primer_Code <- sub("FP_","",ReadData$Primer_Code)
ReadData$Primer_Code <- sub(",[0-9]+","",ReadData$Primer_Code)
ReadData$Primer_Reverse <- sub("RP_","",ReadData$Primer_Reverse)

ReadData$Sample_ID <- paste(substring(ReadData$Acc,1,9),ReadData$Primer_Code,sep="_")
##########################
## End Clipping and preprocessing
##########################

##########################
## Store Data in A SQLite DB
##########################

library(RSQLite)
drv <- SQLite()
con <- dbConnect(drv, dbname="readdata.sqlite")

sql <- "INSERT INTO readdata VALUES ($Acc, $Run, $Sample_ID, $RawLength, $RocheLC, $RocheRC, $RocheLength, 
  $AdapterLC, $AdapterRC, $AdapterLength, $Primer_Code, $Barcode, $FPErr, $Code_Dist, 
  $Primer_Reverse, $RPErr, $keepAdapter, $lucyLC, $lucyRC, $lucyLength, 
  $LucyUnique, $lucyNs, $lucymHomoPrun, $LucyFlip, $LucyRDPgenus, $LucyRDPboot, 
  $LucyTE, $LucyQM, $keepLucy)"


dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = ReadData)

dbGetQuery(con, "SELECT * FROM readdata LIMIT 3")

### Clean up
system(paste("mv TMP.rdp.align ", basefilename, ".lucy.align",sep=""))
system("rm -rf TMP*")
system("rm -rf mothur*")

## THE END :)


