###############################################################################
###
### Periferal function for Amplicon Processing
###
###############################################################################

###### Functions Used During Preprocessing

#########################################################
## check the system for necessary 3rd party applications
#########################################################
## Add check for library files, rdp and mothur go through script and check for others
"checksystem" <- function(){
  print("Checking for required system applications")
  stopifnot(all(
    system("which sfffile") ==0,
    system("which cross_match") ==0,
    system("which lucy") ==0
))}

## trimmed standard deviation
sd.trim <- function(x, trim=0, na.rm=FALSE, ...)
{
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!is.numeric(trim) || length(trim) != 1)
    stop("'trim' must be numeric of length one")
  n <- length(x)
  if(trim > 0 && n > 0) {
    if(is.complex(x)) stop("trimmed sd are not defined for complex data")
    if(trim >= 0.5) return(0)
    lo <- floor(n * trim) + 1
    hi <- n + 1 - lo
    x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}


### read in and format cross_match output
"parse_cm" <- 
  function(filename)
  {
    lines <- readLines(filename)
    align <- grep("^ALIGNMENT",lines,value=T)
    align <- sub(" C ", " ",align)
    #    cm_out <- read.table(filename)
    cm_out <- matrix(unlist(strsplit(align,split=" +")),ncol=13,byrow=T)
    cm_out <- data.frame(cm_out,stringsAsFactors=FALSE)
    cm_out$FC <- "F"
    cm_out$FC[grep("(",cm_out$X11,fixed=TRUE)] <- "C"
    cm_out[cm_out$FC == "C",c("X11","X12","X13")] <- cm_out[cm_out$FC == "C",c("X13","X12","X11")]
    colnames(cm_out) <- c("Alignment","score","perc_sub","perc_del","perc_ins","read_id",
                          "read_start","read_end","read_remain","adapt","adapt_start","adapt_end","adapt_remain","FC")   
    cm_out$read_start <- as.numeric(cm_out$read_start)
    cm_out$read_end <- as.numeric(cm_out$read_end)
    cm_out$read_remain <- as.numeric(gsub("[()]","",cm_out$read_remain))
    cm_out$adapt_remain <- as.numeric(gsub("[()]","",cm_out$adapt_remain))
    cm_out$adapt_start <- as.numeric(cm_out$adapt_start)
    cm_out$adapt_end <- as.numeric(cm_out$adapt_end)
    cm_out$score <- as.numeric(cm_out$score)
    cm_out$perc_sub <- as.numeric(cm_out$perc_sub)
    cm_out$perc_del <- as.numeric(cm_out$perc_del)
    cm_out$perc_ins <- as.numeric(cm_out$perc_ins)
    
    cm_out$perc_sub <- round(cm_out$perc_sub *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
    cm_out$perc_del <- round(cm_out$perc_del *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
    cm_out$perc_ins <- round(cm_out$perc_ins *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
    
    cm_out$read_len <- cm_out$read_end+cm_out$read_remain
    
    cm_out$err <- apply(cm_out[,c("perc_sub","perc_del","perc_ins","adapt_remain","adapt_start")],1,sum) -1
    cm_out
  }


### run the Knitr report generator

##knit(paste("Reports/Preproc_report_",basefilename,".Rmd",sep=""),output=paste("Reports/Preproc_report_",basefilename,".md",sep=""))


### Some figures from the old Preproc_MB.R file

### alignment figures
#pdf(file.path(figurepath,paste(basefilename,".templateEndDist.lucy.pdf",sep="")),width=7.5,height=10)
#hist(align.report$TemplateEnd,breaks=500)
#abline(v=TEtm-4*TEsd,col="blue")
#abline(v=TEtm+4*TEsd,col="blue")
#abline(v=TEtm-75,col="red")
#abline(v=TEtm+75,col="red")
#dev.off()

#pdf(file.path(figurepath,paste(basefilename,".templateEndDistByGENUS.lucy.pdf",sep="")),width=60,height=10)
#boxplot(align.report$TemplateEnd[expand] ~ ReadData$LucyRDPgenus,xaxt="n",xlab="")
#axis(1,labels=FALSE)
#labels <- sort(unique(ReadData$LucyRDPgenus))
#text(1:length(labels), par("usr")[3] - 0.50, srt = 90, adj = 1,
#     labels = labels, xpd = TRUE)
#abline(h=TEtm-75,col="red")
#abline(h=TEtm+75,col="red")
#dev.off()


### Sequence logos
#library(seqLogo)
#lucy <-  read.DNAStringSet(paste(basefilename,".lucy.unique.fasta",sep=""))
#ac <- alphabetByCycle(lucy,alphabet=c("A","C","T","G"))
#acp <- sweep(ac,2,colSums(ac),"/")

#png("logo.lucy.png",height=7,width=30,units="in",pointsize=10,res=300)
#seqLogo(makePWM(acp))
#dev.off()

#adapter <-  read.DNAStringSet(paste(basefilename,".adapter.unique.fasta",sep=""))
#ac <- alphabetByCycle(adapter,alphabet=c("A","C","T","G"))
#acp <- sweep(ac,2,colSums(ac),"/")

#png("logo.adapter.png",height=7,width=30,units="in",pointsize=10,res=300)
#seqLogo(makePWM(acp))
#dev.off()

## Histogram of read cuts
# d <- density(ReadData$RawLength)
# d2 <- density(ReadData$RocheLength)
# d3 <- density(ReadData$AdapterLength[ReadData$keepAdapter])
# d4 <- density(ReadData$lucyLength[ReadData$keepLucy])
# 
# pdf(file.path(figurepath,paste(basefilename,"readlengths.pdf",sep=".")),bg="white")
# plot(d, main="histogram of read lengths")
# lines(d2,col="red")
# lines(d3,col="green")
# lines(d4,col="orange")
# legend("topright",legend=c("Raw","Roche","Primer Trimmed","Lucy Trimmed"),pch=4,text.col=c("black","red","green","orange"),col=c("black","red","green","orange"))
# dev.off()
# 
