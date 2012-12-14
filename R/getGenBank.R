## testing
#access.nb=unique(Lact.align$TemplateName); seq.names = access.nb; species.names = TRUE; 
#gene.names = FALSE; as.character = FALSE

"get.GenBank" <- 
function (access.nb, seq.names = access.nb, species.names = TRUE, 
          definition = TRUE, as.character = FALSE) 
{
  require("Biostrings")
  N <- length(access.nb)
  nrequest <- N%/%400 + as.logical(N%%400)
  X <- character(0)
  for (i in 1:nrequest) {
    a <- (i - 1) * 400 + 1
    b <- 400 * i
    if (i == nrequest) 
      b <- N
    URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                 paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
                 sep = "")
    X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
  }
  VE <- grep("^VERSION", X)
  FI <- grep("^ {0,}ORIGIN", X) + 1
  LA <- which(X == "//") - 1
  obj <- character(length(VE))
  for (i in 1:length(VE)) {
    tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
    obj[i] <- paste(tmp,collapse="")
  }
  GB <- sapply(strsplit(X[VE],split=" +"),"[[",2L)
  GI <- sapply(strsplit(X[VE],split=" +"),"[[",3L)
  
  names(obj) <-GB
  if (!as.character) 
    obj <- DNAStringSet(obj)
  if (species.names) {
    tmp <- character(length(VE))
    sp <- grep("ORGANISM", X)
    for (i in 1:length(VE)) tmp[i] <- unlist(strsplit(X[sp[i]], " +ORGANISM +"))[2]
    attr(obj, "species") <- gsub(" ", "_", tmp)
  }
  if (definition) {
    sp <- grep("^DEFINITION", X,value=T)
    tmp <- sapply(strsplit(sp,split="  "),"[[",2L)
    attr(obj, "definition") <- gsub("\"$", "", tmp)
  }
  obj
}
