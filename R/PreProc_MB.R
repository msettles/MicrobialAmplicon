#!/mnt/home/msettles/opt/bin/Rscript --vanilla 
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

### Microbial Processing Home
microbe.amplicon.home <- "/mnt/home/msettles/CodeProjects/Rpackages/MicrobialAmplicon"
source(file.path(microbe.amplicon.home,"R","functions.R"))
#source(file.path(microbe.amplicon.home,"R","PreProc_functions.R")) ## in update
screenfile <- file.path(microbe.amplicon.home,"ext.data","screen_27f-534r.combined.fa")

sfffiles <- commandArgs(TRUE)
#sfffiles <- "Amplicon_SFFfiles_AmpProc/HVCR0MB02.sff"

tagkey <- "^FP_"
tag_nucs <- 8
primerkey <- "^RP_"

## filter
maxNs <- 2 ## max number of Ns post lucy filtering
maxforwardprimererrors <- 2 ## max number of errors allowed in the forward primer
maxhammingdisttag <- 1 ## max hamming distance allowed for barcode
minlength <- 100  ## minimum length of acceptable sequence
maxlength <- 600  ## maximum length of sequence allowed
maxhomopol <- 10  ## maximum homopolymer allowed

rdpPath <- "/mnt/home/msettles/opt/rdp_classifier_2.5/rdp_classifier-2.5.jar"
nproc=12
## for speciation
reverseSeq <- TRUE
mothur.template="/mnt/home/msettles/projects/Forney/Bacterial_16S/Alignment_db/silva.bacteria.fasta" 


#### This script preprocesses Roche 454 Reads
## START Rdev

basefilename <- sub(".sff","",basename(sfffiles))
checksystem()

fq <- readSff(sfffiles)
clipMode(fq) <- "raw"

#####################
### Base Roche stats and clip
ReadData <- data.frame(Acc=as.character(id(fq)),
                       Run=substring(as.character(id(fq)),1,9),
                       RawLength=width(fq))

ReadData$RocheLC <- start(qualityClip(fq))
ReadData$RocheRC <- end(qualityClip(fq))

ReadData$RocheLength <- width(sread(fq,clipmode="full"))

### End Base Roche
#####################
###
#####################
### Run cross_match looking for primers and tags
#cross_match test.fasta /mnt/home/msettles/CodeProjects/Rpackages/MicrobialAmplicon/ext.data/screen_27f-534r.combined.fa -minmatch 8 -minscore 16 -tags > test.cmout

### make faster by using multicore and splitting reads
writeFastaQual(fq,"TMP.raw",append=FALSE)
system(paste("cross_match TMP.raw.fasta ", screenfile, " -minmatch 8 -minscore 12 -tags > TMP.cmout",sep=""))
cm_out <- parse_cm("TMP.cmout")

cm_out <- cm_out[-intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "C")),]
cm_out <- cm_out[-intersect(grep(primerkey,cm_out$adapt),which(cm_out$FC == "F")),]

### End CrossMatch
############################
###
############################
### Adapter Clip Points
## Ignore Roche's clip points
ReadData$AdapterLC <- 5

### SET new clip points, recording primer information

#### New Left Clip point starts after the Forward Primer, which begins at 1
newLC = cm_out$read_end[(cm_out$FC == "F" & cm_out$read_start == 1)] + cm_out$adapt_remain[(cm_out$FC == "F" & cm_out$read_start == 1)] + 1

tmp <-tapply(newLC,cm_out$read_id[(cm_out$FC == "F" & cm_out$read_start== 1)],max)
ReadData$AdapterLC[match(names(tmp),ReadData$Acc)] <-
  pmax(ReadData$AdapterLC[match(names(tmp),ReadData$Acc)],tmp)

ReadData$AdapterLC <- pmin(ReadData$RawLength,ReadData$AdapterLC)

## Calculate Primer errors INCLUDES TAG!!
primer_offset <-  as.numeric(names(sort(table(cm_out$adapt_start[which(cm_out$FC == "F" & cm_out$read_start == 1)]),decreasing=T))[1])
ReadData$Barcode <- NA
ReadData$Barcode[match(cm_out$read_id[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1))],ReadData$Acc)] <- cm_out$adapt[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1))]
ReadData$FPErr <- 100
ReadData$FPErr[match(cm_out$read_id[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1))],ReadData$Acc)] <- 
  cm_out[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1)),"err"] - primer_offset +1

## Computing, primer tag hamming distance from assigned tag
tags <- subseq(screen,primer_offset+4,width=tag_nucs)[match(ReadData$Barcode[!is.na(ReadData$Barcode)],names(screen))]
rtags <- subseq(sread(fq)[!is.na(ReadData$Barcode)],5,width=tag_nucs)
tagpairs <- cbind(as.character(tags),as.character(rtags))
tagdiff <- numeric(nrow(tagpairs))  ## Sets all to 0
tagdiff[tagpairs[,1] != tagpairs[,2]] <- apply(tagpairs[tagpairs[,1] != tagpairs[,2],],1,stringDist, method="hamming")
ReadData$Code_Dist[!is.na(ReadData$Barcode)] <- tagdiff

##########################################################
#### New Right Clip points occurs before the inner most Reverse Primers
ReadData$AdapterRC <- ReadData$RawLength
## 2 times the remaining adapter to allow for additional homopolymer errors

cm_out2 <- cm_out[order(cm_out$read_id,cm_out$FC,-cm_out$score,cm_out$read_start),]
cm_out2 <- cm_out2[(cm_out2$FC == "C"),]
cm_out2 <- cm_out2[!duplicated(cm_out2$read_id),]
newRC = cm_out2$read_start - 2*cm_out2$adapt_remain - 1

ReadData$AdapterRC[match(cm_out2$read_id,ReadData$Acc)] <-
  pmin(ReadData$AdapterRC[match(cm_out2$read_id,ReadData$Acc)],newRC)

ReadData$AdapterRC <- pmax(5,ReadData$AdapterRC)
ReadData$AdapterLC <- pmin(ReadData$AdapterLC,ReadData$AdapterRC)


# Determining end adapter presence and error
ReadData$Primer_3prime <- NA
ReadData$Primer_3prime[match(cm_out2$read_id[grep(primerkey,cm_out2$adapt)],ReadData$Acc)] <- cm_out2$adapt[grep(primerkey,cm_out2$adapt)]
ReadData$RPErr <- 100
ReadData$RPErr[match(cm_out2$read_id[grep(primerkey,cm_out2$adapt)],ReadData$Acc)] <- 
  cm_out2[grep(primerkey,cm_out2$adapt),"err"]

ReadData$AdapterLength <- ReadData$AdapterRC - ReadData$AdapterLC +1
### End Adapter Assignment
##########################

### Write out data
primerRange <- IRanges(start=ReadData$AdapterLC,ReadData$AdapterRC)
customClip(fq) <- primerRange ## need to add replacement function 
clipMode(fq)  <- "custom"
writeFastaQual(fq,"TMP.adapterClip",append=FALSE)

###########################
## check for poor regions with lucy
# Performing lucy based quality trimming, error (27)
system(paste("lucy -xtra", nproc,"-minimum 0 -debug TMP.lucy_clip.txt -error 0.002 0.002 -output TMP.lucy.fasta TMP.lucy.fasta.qual  TMP.adapterClip.fasta  TMP.adapterClip.fasta.qual"))

lucy <- read.table("TMP.lucy_clip.txt",as.is=T)

ReadData$lucyLC <- ReadData$AdapterLC
ReadData$lucyRC <- ReadData$AdapterRC

zeros <- which(lucy$V9 == 0)
ReadData$lucyRC[match(lucy[zeros,1] ,ReadData$Acc)] <- ReadData$lucyLC[match(lucy[zeros,1] ,ReadData$Acc)]
ReadData$lucyRC[match(lucy[-zeros,1],ReadData$Acc)] <- ReadData$AdapterLC[match(lucy[-zeros,1],ReadData$Acc)] + lucy[-zeros,10] -1
ReadData$lucyLC[match(lucy[-zeros,1],ReadData$Acc)] <- ReadData$AdapterLC[match(lucy[-zeros,1],ReadData$Acc)] + lucy[-zeros,9] -1
ReadData$lucyLength <- ReadData$lucyRC - ReadData$lucyLC +1

# computer read statistics
lucyRange <- IRanges(start=ReadData$lucyLC,ReadData$lucyRC)
customClip <- lucyRange
clipMode(fq)  <- "custom"

fa.lucy <- sread(fq)
## get unique sequences
ReadData$lucyUnique <- NA
fa.lucy <- sort(fa.lucy)
lucy.unique <- length(unique(fa.lucy))
UniqueID <- rep(paste(basefilename,seq.int(lucy.unique),sep="."),times=diff(c(which(!duplicated(fa.lucy)),length(fa.lucy)+1)))
ReadData$lucyUnique[match(names(fa.lucy),ReadData$Acc)] <- UniqueID

expand <- match(ReadData$lucyUnique,unique(UniqueID))

lucyAlphFreq <- alphabetFrequency(fa.lucy[!duplicated(UniqueID)])

ReadData$lucyNs <- lucyAlphFreq[expand,"N"]
homoRuns <-  sapply(fa.lucy[!duplicated(UniqueID)],function(x) sapply(c("A","C","T","G","N"),function(y) longestConsecutive(as.character(x),letter=y)))
ReadData$lucymHomoPrun <- apply(t(homoRuns),1,max)[expand]

######################
## Run the RDP classifier and align using mothur on all the data
######### Lucy Clipped

primerRange <- IRanges(start=ReadData$lucyLC,ReadData$lucyRC)
customClip(fq) <- primerRange ## need to add replacement function 
clipMode(fq)  <- "custom"

fq_uni <- fq[(!duplicated(ReadData$lucyUnique))]
names(fq_uni) <- ReadData$lucyUnique[(!duplicated(ReadData$lucyUnique))]

writeFastaQual(fq_uni,"TMP.rdp",append=FALSE)

######################
system(paste("java -Xmx1g -jar ",rdpPath," -q TMP.rdp.fasta -o TMP.lucy.rdpV6.fix -f fixrank",sep=""))

rdp.lucy <- read.table("TMP.lucy.rdpV6.fix",sep="\t")
rdp.lucy <- rdp.lucy[match(unique(ReadData$lucyUnique),rdp.lucy[,1]),]
rdp.lucy[,1] <- unique(ReadData$lucyUnique)

flip <- read.table("TMP.rdp.flip.accnos",sep="\t")
if(ncol(flip) > 1){
  rdp.lucy$flip <- FALSE
  rdp.lucy$flip[match(flip[grep("reverse complement produced a better alignment",flip[,2]),1],rdp.lucy[,1])] <- TRUE
}

system(paste("mothur \"#align.seqs(candidate=TMP.rdp.fasta, template=", mothur.template ,", flip=T, processors=",nproc,")\"",sep=""))
system(paste("mothur \"#filter.seqs(fasta=TMP.rdp.align, processors=",nproc,");\"",sep=""))

align.report <- read.table("TMP.rdp.align.report",sep="\t",header=T)
align.report <- align.report[match(unique(ReadData$lucyUnique),align.report[,1]),]
align.report[,1] <- unique(ReadData$lucyUnique)

expand <- match(ReadData$lucyUnique,align.report[,1]) ## do I need or use?

TEtm <- 534-27
ReadData$keep <- FALSE
ReadData$keep <-    
  align.report$QueryStart[expand] < 5 &
  align.report$QueryLength[expand]-((align.report$QueryEnd-align.report$QueryStart)[expand]+1) < 5 &
  rdp.lucy$flip[expand] == reverseSeq &

  align.report$TemplateEnd[expand] > (TEtm - 75) & 
  align.report$TemplateEnd[expand] < (TEtm + 75) &

  ReadData$FPErr <= maxforwardprimererrors &
  ReadData$Code_Dist <= maxhammingdisttag &
  ReadData$lucyNs <= maxNs &
  ReadData$lucymHomoPrun <= maxhomopol &
  ReadData$lucyLength > minlength &
  ReadData$lucyLength < maxlength &
  !is.na(ReadData$Barcode)

ReadData$keep[is.na(ReadData$keep)] <- FALSE


### clean up primer designations
ReadData$Barcode <- sub("FP_","",ReadData$Barcode)
ReadData$Barcode <- sub(",[0-9]+","",ReadData$Barcode)
ReadData$Barcode <- sub("FP_","",ReadData$Barcode)
ReadData$Barcode <- sub(",[0-9]+","",ReadData$Barcode)
ReadData$Primer_3prime <- sub("RP_","",ReadData$Primer_3prime)

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
  $AdapterLC, $AdapterRC, $AdapterLength, $Barcode, $FPErr, $Code_Dist, 
  $Primer_3prime, $RPErr, $lucyLC, $lucyRC, $lucyLength, $lucyUnique, $lucyNs, $lucymHomoPrun,
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
###### TO DO: Need to add in Run
sql <- "INSERT INTO rdp_report VALUES ($V1, RUN, $V2, $V3, $V5, $V6, $V8, $V9, $V11, $V12, $V14, $V15, $V17, $V18, $V20, 'NA', 'NA')"

dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = rdp.lucy)
dbCommit(con)

system(paste("mv TMP.rdp.align Output_Files/",basefilename,".lucy.align",sep=""))
system("rm -rf TMP*")
system("rm -rf mothur*")


