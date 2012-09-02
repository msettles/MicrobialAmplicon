#!/mnt/home/msettles/opt/bin/Rscript --vanilla --verbose
###############################################################################
###
### Primary functions/main for Amplicon Processing
###
###############################################################################
options(width=120)
options("stringsAsFactors" = FALSE)
library(rSFFreader)
library(getopt)

version = "Rcode:1.0;rdp:2.5;mothur:1.27;alignment_db:silva.bacteria.fasta"
source("functions.R")

### Microbial Processing Home
microbe.amplicon.home <- "/mnt/home/msettles/CodeProjects/Rpackages/MicrobialAmplicon"
version = "Rcode:1.0;rdp:2.5;mothur:1.27;alignment_db:silva.bacteria.fasta"
source(file.path(microbe.amplicon.home,"R","functions.R"))
screenfile <- file.path(microbe.amplicon.home,"ext.data","screen_27f-534r.combined.fa")

sfffiles <- commandArgs(TRUE)
#sfffiles <- "../Amplicon_SFFfiles_AmpProc/HRLK7U402.sff"

tagkey <- "^FP_"
tag_nucs <- 8
primerkey <- "^RP_"

## filter
maxNs <- 2
maxforwardprimererrors <- 2
maxhammingdisttag <- 1
minlength <- 100
maxlength <- 600
maxhomopol <- 10

#rdpPath <- "/mnt/home/msettles/opt/rdp_classifier_2.2/rdp_classifier-2.2.jar"
rdpPath <- "/mnt/home/msettles/opt/rdp_classifier_2.5/rdp_classifier-2.5.jar"
nproc=12
## for speciation
reverseSeq <- TRUE
mothur.template="/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta" 


#### This script preprocesses Roche 454 Reads
## START Rdev

basefilename <- sub(".sff","",basename(sfffiles))
checksystem()

outfile <- file.path("Output_Files",paste(basefilename,".out",sep=""))


fq <- readsff(sfffiles)
cat(paste("Total Number of Reads in SFF files:",length(fq),"\n"),file=outfile,append=F)    
clipMode(fq) <- "Raw"

#####################
### Base Roche stats and clip
print("Identifying Roche Clip points")
ReadData <- data.frame(Acc=as.character(id(fq)),
                       Run=substring(as.character(id(fq)),1,9),
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

writeFastaQual(fq,"TMP.raw",append=FALSE)
system(paste("cross_match TMP.raw.fasta ", screenfile, " -minmatch 8 -minscore 16 -tags > TMP.cmout",sep=""))
cm_out <- parse_cm("TMP.cmout")

### keep those forward match > 16 score
cm_out <- cm_out[-which(cm_out$score <= 16),]

### End CrossMatch
############################
###
############################
### Adapter Clip Points
## Ignore Roche's clip points
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
### ADD IN DUMMY BARCODE FOR LATER USE
ReadData$FPErr <- 100
ReadData$FPErr[match(cm_out$read_id[grep(tagkey,cm_out$adapt)],ReadData$Acc)] <- 
  cm_out[grep(tagkey,cm_out$adapt),"err"] - primer_offset +1

cat("Forward Primer Errors:\n",file=outfile,append=T)
sink(file=outfile,append=TRUE)
stem(ReadData$FPErr,scale=0.5,width=10)
sink()
cat("\n",file=outfile,append=T)

ReadData$Barcode <- ReadData$Primer_Code


print("Computing, primer tag hamming distance from assigned tag")
screen <- readDNAStringSet(file.path(screenfile))
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
cat(paste("Reads Greater than Min Length > :",minlength, " :: Reads meeting that criteria ::",table(ReadData$AdapterLength > minlength)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Less than Max Length :",maxlength, " :: Reads meeting that criteria ::",table(ReadData$AdapterLength < maxlength)[2], "\n"),file=outfile,append=T)
cat(paste("Primer with Tab found :",table(!is.na(ReadData$Primer_Code))[2], "\n"),file=outfile,append=T)

### Write out data
primerRange <- IRanges(start=ReadData$AdapterLC,ReadData$AdapterRC)
adapterClip(fq) <- primerRange ## need to add replacement function 
fq@adapterClip <- primerRange
clipMode(fq)  <- "Full"

writeFastaQual(fq,"TMP.adapterClip",append=FALSE)

###########################
## check for poor regions with lucy
print("Performing lucy based quality trimming, error (27)")
### also limits size to 100
system(paste("lucy -xtra", nproc,"-minimum 0 -debug TMP.lucy_clip.txt -error 0.002 0.002 -output TMP.lucy.fasta TMP.lucy.fasta.qual  TMP.adapterClip.fasta  TMP.adapterClip.fasta.qual"))

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
ReadData$lucyUnique <- NA
fa.lucy <- sort(fa.lucy)
lucy.unique <- length(unique(fa.lucy))
UniqueID <- rep(paste(basefilename,seq.int(lucy.unique),sep="."),times=diff(c(which(!duplicated(fa.lucy)),length(fa.lucy)+1)))
ReadData$lucyUnique[match(names(fa.lucy),ReadData$Acc)] <- UniqueID

cat(paste("total number of unique reads:", lucy.unique,"\n"),file=outfile,append=T)

expand <- match(ReadData$lucyUnique,unique(UniqueID))

lucyAlphFreq <- alphabetFrequency(fa.lucy[!duplicated(UniqueID)])

ReadData$lucyNs <- lucyAlphFreq[expand,"N"]
homoRuns <-  sapply(fa.lucy[!duplicated(UniqueID)],function(x) sapply(c("A","C","T","G","N"),function(y) longestConsecutive(as.character(x),letter=y)))
ReadData$lucymHomoPrun <- apply(t(homoRuns),1,max)[expand]

cat(paste("\nLucy Additional Clipping\n"),file=outfile,append=T)
cat(paste("Maximum Number of Abiguous Characters (Ns) <= :",maxNs, " :: Reads meeting that criteria ::",table(ReadData$lucyNs <= maxNs)[2], "\n"),file=outfile,append=T)
cat(paste("Maximum Length of a Homopolymer Run <= :",maxhomopol, " :: Reads meeting that criteria ::",table(ReadData$lucymHomoPrun <= maxhomopol)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Greater than Min Length > :",minlength, " :: Reads meeting that criteria ::",table(ReadData$lucyLength > minlength)[2], "\n"),file=outfile,append=T)
cat(paste("Reads Less than Max Length :",maxlength, " :: Reads meeting that criteria ::",table(ReadData$lucyLength < maxlength)[2], "\n"),file=outfile,append=T)


######################
## Run the RDP classifier on all the data
######### Lucy Clipped

fq_uni <- fq[(!duplicated(ReadData$lucyUnique))]
names(fq_uni) <- ReadData$lucyUnique[(!duplicated(ReadData$lucyUnique))]

writeFastaQual(fq_uni,"TMP.rdp",append=FALSE)

######################
system(paste("java -Xmx1g -jar ",rdpPath," -q TMP.rdp.fasta -o TMP.lucy.rdpV6.fix -f fixrank",sep=""))
rdp.lucy <- read.table("TMP.lucy.rdpV6.fix",sep="\t")
rdp.lucy <- rdp.lucy[match(unique(ReadData$lucyUnique),rdp.lucy[,1]),]
rdp.lucy[,1] <- unique(ReadData$lucyUnique)

system(paste("mothur \"#align.seqs(candidate=TMP.rdp.fasta, template=", mothur.template ,", flip=T, processors=",nproc,")\"",sep=""))

flip <- read.table("TMP.rdp.flip.accnos",sep="\t")
align.report <- read.table("TMP.rdp.align.report",sep="\t",header=T)
align.report <- align.report[match(unique(ReadData$lucyUnique),align.report[,1]),]
align.report[,1] <- unique(ReadData$lucyUnique)
if(ncol(flip) > 1){
  align.report$flip <- FALSE
  align.report$flip[match(flip[grep("reverse complement produced a better alignment",flip[,2]),1],align.report[,1])] <- TRUE
}
align.report$QueryFull <- ReadData$lucyLength[match(align.report[,1],ReadData$lucyUnique)]
align.report$adpLC <- ReadData$AdapterLC[match(align.report[,1],ReadData$lucyUnique)]

expand <- match(ReadData$lucyUnique,align.report[,1]) ## do I need or use?

TEtm <- 534-27
ReadData$keep <- FALSE
ReadData$keep <-    
  align.report$QueryStart[expand] < 5 &
  (align.report$QueryFull-align.report$QueryLength)[expand] < 5 &
  align.report$flip[expand] == reverseSeq &
  align.report$TemplateEnd[expand] > (TEtm - 75) & 
  align.report$TemplateEnd[expand] < (TEtm + 75) &
  ReadData$FPErr <= maxforwardprimererrors &
  ReadData$Code_Dist <= maxhammingdisttag &
  ReadData$lucyNs <= maxNs &
  ReadData$lucymHomoPrun <= maxhomopol &
  ReadData$lucyLength > minlength &
  ReadData$lucyLength < maxlength &
  !is.na(ReadData$Primer_Code)

ReadData$keep[is.na(ReadData$keep)] <- FALSE

cat(paste("RESULTS\nAdditional Lucy clipped filter: ",sum(ReadData$keep), "Reads kept ::", sum(!ReadData$keep), "Reads removed\n"),file=outfile,append=T)

###########################
cat("Overlap between Adapter Clipped and Lucy Clipped Data:\n",file=outfile,append=T)
cat("\n",file=outfile,append=T)
cat(paste("Mean number of bases removed by Lucy clipping:",signif(mean(ReadData$AdapterLength-ReadData$lucyLength),3),"\n"),file=outfile,append=T)
cat(paste("Median number of bases removed by Lucy clipping:",signif(median(ReadData$AdapterLength-ReadData$lucyLength),3),"\n"),file=outfile,append=T)
cat(paste("Number of reads Lucy    Clipped and Unique:",table(!duplicated(ReadData$lucyUnique) & ReadData$keep)[2],"\n"),file=outfile,append=T)
##########################

### clean up primer designations
ReadData$Primer_Code <- sub("FP_","",ReadData$Primer_Code)
ReadData$Primer_Code <- sub(",[0-9]+","",ReadData$Primer_Code)
ReadData$Barcode <- sub("FP_","",ReadData$Barcode)
ReadData$Barcode <- sub(",[0-9]+","",ReadData$Barcode)
ReadData$Primer_Reverse <- sub("RP_","",ReadData$Primer_Reverse)

ReadData$version <- version
##########################
## End Clipping and preprocessing
##########################

############################
## Store Data in A SQLite DB
############################

library(RSQLite)
drv <- SQLite()
con <- dbConnect(drv, dbname="amplicondata.sqlite")

sql <- "INSERT INTO read_data VALUES ($Acc, $Run, $RawLength, $RocheLC, $RocheRC, $RocheLength, 
  $AdapterLC, $AdapterRC, $AdapterLength, $Primer_Code, $FPErr, $Barcode, $Code_Dist, 
  $Primer_Reverse, $RPErr, $lucyLC, $lucyRC, $lucyLength, $lucyUnique, $lucyNs, $lucymHomoPrun,
  $keep, $version)"

dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = ReadData)
dbCommit(con)

align.report$SimBtwnQuery <- align.report$SimBtwnQuery.Template
sql <- "INSERT INTO align_report VALUES ($QueryName, $QueryLength, $TemplateName, $TemplateLength, $SearchMethod,
  $SearchScore, $AlignmentMethod, $QueryStart, $QueryEnd, $TemplateStart, $TemplateEnd, $PairwiseAlignmentLength, 
  $GapsInQuery, $GapsInTemplate, $LongestInsert, $SimBtwnQuery, $flip, $QueryFull, $adpLC)"

dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = align.report)
dbCommit(con)

sql <- "INSERT INTO rdp_report VALUES ($V1, $V2, $V3, $V5, $V6, $V8, $V9, $V11, $V12, $V14, $V15, $V17, $V18, $V20, 'NA', 'NA')"


dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = rdp.lucy)
dbCommit(con)

system(paste("mv TMP.rdp.align Output_Files/",basefilename,".lucy.align",sep=""))
system("rm -rf TMP*")
system("rm -rf mothur*")

