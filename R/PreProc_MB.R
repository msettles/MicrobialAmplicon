#!/mnt/home/msettles/opt/bin/Rscript --vanilla 
###############################################################################
###
### Primary functions/main for Amplicon Processing
###
###############################################################################
## setup
options(width=120)
options("stringsAsFactors" = FALSE)
library(RSQLite)
library(rSFFreader)
library(getopt)
library(multicore)
library(knitr)


### Microbial Processing Home
microbe.amplicon.home <- "/mnt/home/msettles/CodeProjects/Rpackages/MicrobialAmplicon"
alignment_db <- "/mnt/home/msettles/Projects/AmpliconProcessing/"
nproc=12
dbname <- "amplicondataV2.0.sqlite"

source(file.path(microbe.amplicon.home,"R","PreProc_functions.R"))

checksystem()


sfffiles <- commandArgs(TRUE)
#sfffiles <- "Amplicon_SFFfiles_AmpProc/HRHT7FQ02.sff"
basefilename <- sub(".sff","",basename(sfffiles))

source(file.path(microbe.amplicon.home,"R","DButilities.R"))
con <- dbCon(dbname=dbname)   ### need to add check for connection
if (any(!is.na(match(basefilename,getProcessedRuns(con))))) stop("SFFfile(s) already processed")

source(file.path(microbe.amplicon.home,"R","PreProc_options.R"))

###########################################################################################
# Adapters, Barcodes and Primers
###########################################################################################

#####################
### Record Base Roche stats and clip
fq <- readSff(sfffiles)
clipMode(fq) <- "raw"

ReadData <- data.frame(Acc=as.character(id(fq)),
                       Run=substring(as.character(id(fq)),1,9),
                       RawLength=width(fq))

ReadData$RocheLC <- start(qualityClip(fq))
ReadData$RocheRC <- end(qualityClip(fq))
ReadData$RocheLength <- width(sread(fq,clipmode="full"))

#####################
### Run cross_match looking for primers and tags

### can make faster by using multicore and splitting reads

### multicore version of crossmatch
mclapply(chunk(seq.int(1,length(fq)),nproc), function(x){
  writeFastaQual(fq[x],paste(basefilename,".",x[1],".raw",sep=""),append=FALSE)
  system(cross_match_call(x[1]))
})

#writeFastaQual(fq,paste(basefilename,".raw",sep=""),append=FALSE)
#system(cross_match_call)

cm_out <- parse_cm(dir(pattern=paste(basefilename,".*.cmout",sep="")))

### Checkpoint
save.image(paste(basefilename,".RData",sep=""))

### remove all impossible match (ie forward primer as compliment)
cm_out <- cm_out[-intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "C")),]
cm_out <- cm_out[-intersect(grep(primerkey,cm_out$adapt),which(cm_out$FC == "F")),]
cm_out <- cm_out[-which(cm_out$FC == "F" & cm_out$score < 20), ] ## forward match must be at least 20
############################
### Adapter Left Clip Points
ReadData$AdapterLC <- 5  ### default start post key

### SET new clip points, recording primer information

#### New Left Clip point starts after the Forward Primer, which begins at 1
newLC = cm_out$read_end[(cm_out$FC == "F" & cm_out$read_start == 1)] + cm_out$adapt_remain[(cm_out$FC == "F" & cm_out$read_start == 1)] + 1

tmp <-tapply(newLC,cm_out$read_id[(cm_out$FC == "F" & cm_out$read_start== 1)],max)
ReadData$AdapterLC[match(names(tmp),ReadData$Acc)] <-
  pmax(ReadData$AdapterLC[match(names(tmp),ReadData$Acc)],tmp)

ReadData$AdapterLC <- pmin(ReadData$RawLength,ReadData$AdapterLC)

## Calculate Primer errors INCLUDES TAG!!
ReadData$Barcode <- NA
ReadData$Barcode[match(cm_out$read_id[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1))],ReadData$Acc)] <- cm_out$adapt[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1))]
#### CHECK THIS CODE "FP_Bender633_CGTTAACG_p534R,8" score of 14

## Computing, primer tag hamming distance from assigned tag
primer_offset <-  as.numeric(names(sort(table(cm_out$adapt_start[which(cm_out$FC == "F" & cm_out$read_start == 1)]),decreasing=T))[1])
screen <- readDNAStringSet(cross_match_screenfile)
tags <- subseq(screen,primer_offset+4,width=tag_nucs)[match(ReadData$Barcode[!is.na(ReadData$Barcode)],names(screen))]
rtags <- subseq(sread(fq)[!is.na(ReadData$Barcode)],5,width=tag_nucs)
tagpairs <- cbind(as.character(tags),as.character(rtags))
tagdiff <- numeric(nrow(tagpairs))  ## Sets all to 0
tagdiff[tagpairs[,1] != tagpairs[,2]] <- apply(tagpairs[tagpairs[,1] != tagpairs[,2],],1,stringDist, method="hamming")
ReadData$Barcode_Err <- NA
ReadData$Barcode_Err[!is.na(ReadData$Barcode)] <- tagdiff
ReadData$FP_Err <- NA
ReadData$FP_Err[match(cm_out$read_id[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1))],ReadData$Acc)] <- 
  cm_out[intersect(grep(tagkey,cm_out$adapt),which(cm_out$FC == "F" & cm_out$read_start == 1)),"err"] - primer_offset +1


############################
### Adapter Right Clip Points
ReadData$AdapterRC <- ReadData$RawLength

### choose match by best score
cm_out2 <- cm_out[(cm_out$FC == "C"),]
cm_out2 <- cm_out2[order(cm_out2$read_id,cm_out2$FC,-cm_out2$score,cm_out2$read_start),]
cm_out2 <- cm_out2[!duplicated(cm_out2$read_id),]
## 2 times the remaining adapter to allow for additional homopolymer errors
newRC = cm_out2$read_start - 2*cm_out2$adapt_remain - 1

ReadData$AdapterRC[match(cm_out2$read_id,ReadData$Acc)] <-
  pmin(ReadData$AdapterRC[match(cm_out2$read_id,ReadData$Acc)],newRC)

ReadData$AdapterRC <- pmax(5,ReadData$AdapterRC)
ReadData$AdapterLC <- pmin(ReadData$AdapterLC,ReadData$AdapterRC)

# Determining end adapter presence and error
ReadData$Primer_3prime <- NA
ReadData$Primer_3prime[match(cm_out2$read_id[grep(primerkey,cm_out2$adapt)],ReadData$Acc)] <- cm_out2$adapt[grep(primerkey,cm_out2$adapt)]
ReadData$RP_Err <- NA
ReadData$RP_Err[match(cm_out2$read_id[grep(primerkey,cm_out2$adapt)],ReadData$Acc)] <- 
  cm_out2[grep(primerkey,cm_out2$adapt),"err"]

ReadData$AdapterLength <- ReadData$AdapterRC - ReadData$AdapterLC +1

############################
### Performing lucy based quality trimming, error (27)

### Write out data
primerRange <- IRanges(start=ReadData$AdapterLC,ReadData$AdapterRC)
customClip(fq) <- primerRange ## need to add replacement function 
clipMode(fq)  <- "custom"
writeFastaQual(fq,paste(basefilename,".adapterClip",sep=""),append=FALSE)

system(lucy_call)

### Checkpoint
save.image(paste(basefilename,".RData",sep=""))

lucy <- read.table(paste(basefilename,".lucy_clip.txt",sep=""),as.is=T)

ReadData$LucyLC <- ReadData$AdapterLC
ReadData$LucyRC <- ReadData$AdapterRC

zeros <- which(lucy$V9 == 0)
ReadData$LucyRC[match(lucy[zeros,1] ,ReadData$Acc)] <- ReadData$LucyLC[match(lucy[zeros,1] ,ReadData$Acc)]
ReadData$LucyRC[match(lucy[-zeros,1],ReadData$Acc)] <- ReadData$AdapterLC[match(lucy[-zeros,1],ReadData$Acc)] + lucy[-zeros,10] -1
ReadData$LucyLC[match(lucy[-zeros,1],ReadData$Acc)] <- ReadData$AdapterLC[match(lucy[-zeros,1],ReadData$Acc)] + lucy[-zeros,9] -1
ReadData$LucyLength <- ReadData$LucyRC - ReadData$LucyLC +1

# compute read statistics
lucyRange <- IRanges(start=ReadData$LucyLC,ReadData$LucyRC)
customClip <- lucyRange
clipMode(fq)  <- "custom"

fa.lucy <- sread(fq)
## get unique sequences
ReadData$LucyUnique <- NA
fa.lucy <- sort(fa.lucy)
lucy.unique <- length(unique(fa.lucy))
UniqueID <- rep(paste(basefilename,seq.int(lucy.unique),sep="."),times=diff(c(which(!duplicated(fa.lucy)),length(fa.lucy)+1)))

ReadData$LucyUnique[match(names(fa.lucy),ReadData$Acc)] <- UniqueID

expand <- match(ReadData$LucyUnique,UniqueID[!duplicated(UniqueID)])

lucyAlphFreq <- alphabetFrequency(fa.lucy[!duplicated(UniqueID)])

ReadData$LucyNs <- lucyAlphFreq[expand,"N"]
## can parallelize to multicore
homoRuns <-  sapply(fa.lucy[!duplicated(UniqueID)],function(x) sapply(c("A","C","T","G","N"),function(y) longestConsecutive(as.character(x),letter=y)))
ReadData$LucymHomoPrun <- apply(t(homoRuns),1,max)[expand]

fq_uni <- fq[(!duplicated(ReadData$LucyUnique))]
names(fq_uni) <- ReadData$LucyUnique[(!duplicated(ReadData$LucyUnique))]

######################
## Run the RDP classifier and align using mothur on all the data
#rdp_call <- function(x) paste("java -Xmx1g -jar ",rdp_path," -q ",basefilename,".",x,".rdp.fasta -o ",basefilename,".",x,".lucy.rdpV6.fix -f fixrank",sep="")

### multicore version of crossmatch

mclapply(chunks <- chunk(seq.int(1,length(fq_uni)),nproc), function(x){
  writeFastaQual(fq_uni[x],paste(basefilename,".",x[1],".rdp",sep=""),append=FALSE)
  system(rdp_call(x[1]))
})

#writeFastaQual(fq_uni,paste(basefilename,".rdp",sep=""),append=FALSE)
#system(rdp_call)


### Checkpoint
save.image(paste(basefilename,".RData",sep=""))

rdp.lucy <- do.call("rbind", lapply(dir(pattern=paste(basefilename,".*.lucy.rdpV6.fix",sep="")), function(fn) 
  read.table(fn,sep="\t")))
                

rdp.lucy$flip <- TRUE
if (any(file.exists(dir(pattern=paste(basefilename,".*.lucy.rdpV6.fix",sep=""))))){
  flip <- do.call("rbind",lapply(dir(pattern=paste(basefilename,".*.lucy.rdpV6.fix",sep="")),function(fn) read.table(fn,sep="\t")))
  if(ncol(flip) > 1){
    rdp.lucy$flip[-match(flip[grep("reverse complement produced a better alignment",flip[,2]),1],rdp.lucy[,1])] <- FALSE
  }
}

## Fills in unique seqs that fail rdp with NAs
rdp.lucy <- rdp.lucy[match(unique(ReadData$LucyUnique),rdp.lucy[,1]),]
rdp.lucy[,1] <- unique(ReadData$LucyUnique)

######################
## Align to reference DB using Mothur
writeFastaQual(fq_uni,paste(basefilename,".rdp",sep=""),append=FALSE)
system(mothur_align_call)
#system(mothur_filter_call)

### Checkpoint
save.image(paste(basefilename,".RData",sep=""))

## Fills in unique seqs that fail alignment with NAs
align.report <- read.table(paste(basefilename,".rdp.align.report",sep=""),sep="\t",header=T)
align.report <- align.report[match(unique(ReadData$LucyUnique),align.report[,1]),]
align.report[,1] <- unique(ReadData$LucyUnique)

######################
## Determine Filters
expand <- match(ReadData$LucyUnique,align.report[,1]) ## do I need or use?

qual_FP_err <- ReadData$FP_Err <= maxforwardprimererrors  & !is.na(ReadData$FP_Err)
qual_hamm_dist <- ReadData$Barcode_Err <= maxhammingdisttag & !is.na(ReadData$Barcode_Err)
qual_maxN <- ReadData$LucyNs <= maxNs
qual_homoPrun <- ReadData$LucymHomoPrun <= maxhomopol

align_start_error <- align.report$QueryStart[expand] <= align_length_max_error
align_length_error <- align.report$QueryLength[expand]-((align.report$QueryEnd-align.report$QueryStart)[expand]+1) <= align_length_max_error
#align_end_diff <- align.report$TemplateEnd[expand] > (TE_exp - TE_max_dist) & 
#  align.report$TemplateEnd[expand] < (TE_exp + TE_max_dist)
rdp_flip <- rdp.lucy$flip[expand] == reverseSeq & !is.na(rdp.lucy$flip[expand])

len_min <- ReadData$LucyLength > minlength
len_max <- ReadData$LucyLength < maxlength


ReadData$keep <- FALSE
ReadData$keep  <-   
  qual_hamm_dist &
  qual_FP_err  &
  qual_maxN &
  qual_homoPrun  &

  align_start_error &
  align_length_error &
#  align_end_diff  &
  rdp_flip  &

  len_min &
  len_max 

ReadData$keep[is.na(ReadData$keep)] <- FALSE

### clean up primer designations
ReadData$Barcode <- sub("FP_","",ReadData$Barcode)
ReadData$Barcode <- sub(",[0-9]+","",ReadData$Barcode)
ReadData$Primer_3prime <- sub("RP_","",ReadData$Primer_3prime)

ReadData$version <- version

### Checkpoint
save.image(paste(basefilename,".RData",sep=""))
############################
## End Clipping and preprocessing
############################

############################
## Generate Knitr Report
############################
setwd("Reports")
file.copy(file.path(microbe.amplicon.home,"report_templates","Preproc_report.Rmd"),file.path(paste("Preproc_report_",basefilename,".Rmd",sep="")),overwrite=TRUE)
knit2html(paste("Preproc_report_",basefilename,".Rmd",sep=""))
setwd("..")
############################
## Store Data in A SQLite DB
############################

dbDis(con)
con <- dbCon(dbname=dbname)

sql <- "INSERT INTO read_data VALUES ($Acc, $Run, $RawLength, $RocheLC, $RocheRC, $RocheLength, 
  $AdapterLC, $Barcode, $Barcode_Err, $FP_Err, $AdapterRC, $Primer_3prime, $RP_Err, $AdapterLength, 
  $LucyLC, $LucyRC, $LucyLength, $LucyUnique, $LucyNs, $LucymHomoPrun,
  $keep, $version)"

dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = ReadData)
dbCommit(con)

align.report$Run <- ReadData$Run[match(align.report$QueryName,ReadData$LucyUnique)]
align.report$SimBtwnQuery <- align.report$SimBtwnQuery.Template
sql <- "INSERT INTO align_report VALUES ($QueryName, $Run, $QueryLength, $TemplateName, $TemplateLength, $SearchMethod,
  $SearchScore, $AlignmentMethod, $QueryStart, $QueryEnd, $TemplateStart, $TemplateEnd, $PairwiseAlignmentLength, 
  $GapsInQuery, $GapsInTemplate, $LongestInsert, $SimBtwnQuery)"

dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = align.report)
dbCommit(con)

rdp.lucy$Run <- ReadData$Run[match(rdp.lucy$V1,ReadData$LucyUnique)]
sql <- "INSERT INTO rdp_report VALUES ($V1, $Run, $flip, $V3, $V5, $V6, $V8, $V9, $V11, $V12, $V14, $V15, $V17, $V18, $V20, 'NA')"

dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = rdp.lucy)
dbCommit(con)

### Update Completed Run table
updateProcessedRuns(con,unique(ReadData$Run),sum(ReadData$keep),nrow(ReadData)-sum(ReadData$keep))

dbDis(con)
#### cleanup
system(paste("rm -rf ",basefilename,"*",sep=""))
system("rm -rf mothur*")
system("rm -rf .RData .Rhistory")

