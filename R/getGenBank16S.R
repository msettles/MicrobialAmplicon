## testing
taxonomicID="Lactobacillaceae"
searchTerm="AND (16S ribosomal RNA) NOT (Uncultured)"
"get.GenBank16s" <- 
function (taxonomicID, searchTerm) 
{
  require("Biostrings")
  require("stringr")
  N <- length(names)
  X <- vector("list",N)
  for (i in 1:length(N)) {
    URL <- paste("http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&retmax=1&usehistory=y&term=(",taxonomicID,"%5BOrganism%5D) ",searchTerm,sep="")
    URL <- gsub(" ","%20",URL)
    Y <- scan(file=URL, what="",sep="\n",quiet=TRUE)
    named.p <-
      paste("m|<Count>(\\d+)</Count>.*<QueryKey>(\\d+)</QueryKey>.*<WebEnv>(\\S+)</WebEnv>|s;",sep="")
    keys <- str_match_all(Y[3],named.p)[[1]]
    Counts <- as.numeric(keys[1,2])
    Query_Key <- keys[1,3]
    WebEnv <- keys[1,4]
    print(paste("Retrieving ", Counts," sequences",sep=""))
    retstart = 0
    retmax=400
    rettype="gb"
    db = "nucleotide"
    nrequest <- Counts%/%retmax + as.logical(Counts%%retmax)
    Z <- character(0)
    for (j in 1:nrequest) {
#    for (j in 1:10) {
        a <- (j - 1) * retmax
      b <- retmax * j
      if (j == nrequest) 
        b <- Counts
      URL <- paste("http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                      "rettype=",rettype,"&retmode=text&retstart=",a,"&retmax=",retmax,"&",
                      "db=",db,"&query_key=",Query_Key,"&WebEnv=",WebEnv,sep="")
      Z <- c(Z, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
      cat("."); flush.console()
    }
    cat("\n"); flush.console()
    X[[i]] <- Z
  }
  X <- X[[1]] ### DOES NOT WORK WITH MULTIPLE TERMS YET

  is.na(X) <- X == "//"
  R1 <- rle(!is.na(X))
  X <- split(X, rep(cumsum(R1$values) * R1$values, R1$lengths))[-1]
  X <- X[as.logical(sapply(X,function(x) length(grep("^ {0,}ORIGIN", x))))] ## REMOVE THOSE WITH NO SEQUENCE
  
  VERSION <- sapply(X, function(x) {
    VE <- grep("^VERSION", x,value=T)
    GB <- sapply(strsplit(VE,split=" +"),"[[",2L)
#    GI <- sapply(strsplit(VE,split=" +"),"[[",3L)
    GB
  })

  SEQUENCE <- sapply(X, function(x){
    FI <- grep("^ {0,}ORIGIN", x) + 1    
    tmp <- gsub("[[:digit:] ]", "", x[FI:length(x)])
    obj <- paste(tmp,collapse="")
    #  grep("gene|tRNA|rRNA|CDS",t)
    #  grep("product",t)
    obj
  })
  names(SEQUENCE) <-VERSION
  SEQUENCE <- DNAStringSet(SEQUENCE)

  ORGANISM <- sapply(X, function(x){
    ORG <- grep("^ {0,}ORGANISM", x)
    END <- grep("^REFERENCE",x)[1] -1
    tmp <- gsub("[ .]", "", x[(ORG+1):END])
    obj <- paste(tmp,collapse="")
    strain <- gsub("^[ ]+ORGANISM +","",x[ORG])
    species <- paste(strsplit(strain,split=" ")[[1]][1],strsplit(strain,split=" ")[[1]][2])
    obj <- paste(obj,species,strain,sep=";")
  })
  
  DEFINITION <- sapply(X, function(x){
    sp <- grep("^DEFINITION", x,value=T)
    tmp <- sapply(strsplit(sp,split="  "),"[[",2L)
    tmp
  })

  attr(SEQUENCE,"taxonomy") <- ORGANISM
  attr(SEQUENCE,"description") <- DEFINITION
  SEQUENCE
}
