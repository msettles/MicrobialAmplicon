#####################################################
####  Functions based on DB
#####################################################


"dbCon" <- function(dbname="amplicondata.sqlite"){
  require(RSQLite)
  drv <- SQLite()
  dbConnect(drv, dbname=dbname)
}

"dbDis" <- function(con) dbDisconnect(con)


"getProcessedRuns" <- function(con){
  dbGetQuery(con, "SELECT DISTINCT Run FROM read_data;")$Run
}

"getProjects" <- function(con){
  dbGetQuery(con, "SELECT DISTINCT Project FROM pool_metadata;")$Project
}

"getRuns" <- function(con){
  dbGetQuery(con, "SELECT DISTINCT Run FROM read_data;")$Run
}

## pool_mapping to add
"updateRunMappings" <- function(con,filename="MetaData/Pool_Mapping.txt"){
  ## get the current table
  mappings <- dbReadTable(con,"pool_mapping")
  pooltable <- read.table(filename,sep="\t")
  pooltable <- pooltable[!(apply(pooltable,1,paste,collapse=" ") %in% apply(mappings,1,paste,collapse=" ")),]
  print(paste("Adding ",nrow(pooltable)," entries",sep=""))
  if (nrow(pooltable) > 0)
    dbWriteTable(con,"pool_mapping",pooltable,row.names=F,append=T)
}

## Metadata Table to add
"updateMetaData" <- function(con,filename="MetaData/Pool_MetaData_Sept28-2012.txt"){
  # get the current table
  metadata <- dbReadTable(con,"pool_metadata")
  metatable <- read.table(filename,sep="\t",header=T)
  metatable <- metatable[!(apply(metatable[,c("Project","Sample_ID","Reverse_Primer")],1,paste,collapse=" ") %in% apply(metadata[,c("Project","Sample_ID","Reverse_Primer")],1,paste,collapse=" ")),]
  print(paste("Adding ",nrow(metatable)," entries",sep=""))
  if (nrow(pooltable) > 0)
    dbWriteTable(con,"pool_metadata",metatable,row.names=F,append=T)
}

"sample_read_counts" <- function(con,project,outfile, dir="OutputFiles"){

  if (missing(outfile)) outfile <- paste(project,".readcounts.txt",sep="")

  reads <- dbGetQuery(con,
                         paste("Select pool_metadata.Sample_ID, read_data.* 
                                FROM pool_metadata, read_data, pool_mapping 
                                WHERE pool_metadata.project='",project,"' 
                                  AND pool_metadata.Pool=pool_mapping.Pool 
                                  AND pool_metadata.Reverse_Primer=read_data.Primer_Code 
                                  AND pool_mapping.Run=read_data.Run",sep=""))
  samples <-dbGetQuery(con,
                         paste("Select Run, pool_metadata.Pool, Reverse_Primer, Sample_ID, ID 
                                FROM pool_metadata, pool_mapping 
                                WHERE pool_metadata.Pool=pool_mapping.Pool 
                                  AND pool_metadata.project='",project,"'",sep=""))
  read_counts <- table(paste(reads$Run,reads$Barcode,reads$Sample_ID),reads$keep)
  colnames(read_counts) <- c("Fail","Pass")
  rcTable <- data.frame(samples,read_counts[match(paste(samples$Run,samples$Reverse_Primer,samples$Sample_ID),rownames(read_counts)),])
  
  dir.create(path=file.path(dir,project),recursive=T,showWarnings=F)
  write.table(rcTable,file=file.path(dir,project,outfile),sep="\t",quote=F,col.names=T,row.names=F)

}

### computes abundance tables, based on project
"abundance.table" <- function(con,project,outfile, dir="OutputFiles", rdpThres = 0.5){

  if (missing(outfile)) outfile <- paste(project,".abundance.txt",sep="")

  samples <-dbGetQuery(con,
                         paste("Select Sample_ID
                                FROM pool_metadata, pool_mapping 
                                WHERE pool_metadata.Pool=pool_mapping.Pool 
                                  AND pool_metadata.project='",project,"'",sep=""))

  rdp.otu <- dbGetQuery(con,
                         paste("Select pool_metadata.Sample_ID, rdp_report.*, read_data.keep 
                                FROM pool_metadata, read_data, pool_mapping, rdp_report 
                                WHERE pool_metadata.Project='",project,"' 
                                  AND pool_metadata.Pool=pool_mapping.Pool 
                                  AND pool_metadata.Reverse_Primer=read_data.Primer_Code 
                                  AND pool_mapping.Run=read_data.Run 
                                  AND read_data.lucyUnique=rdp_report.QueryName",sep=""))

  rdp.otu <- rdp.otu[rdp.otu$keep == 1,]

  numberoflevels <- (ncol(rdp.otu)-4)/2
  abtable <- data.frame(id=rdp.otu[,"Sample_ID",],assign=as.character(rdp.otu$domain_name),level="domain")
  
  taxa_levels <- c("domain","phylum","class","order","family","genus","species")
  abtable$assign[(rdp.otu$domain_bootstrap < rdpThres)] <- "Unclassified"
  abtable$level[(rdp.otu$domain_bootstrap < rdpThres)] <- "Unknown"
  
  for (i in seq.int(2,numberoflevels)){
    column = 2*i+2
    prows <- rdp.otu[, (column + 1)] >= rdpThres & (rdp.otu[, (column + 1)] != "NA")
    
    abtable$assign[prows] <- as.character(rdp.otu[,column][prows])
    abtable$level[prows] <- taxa_levels[i]
  }
  
  
  abundanceTable <- table(abtable$assign,abtable$id)
  abundanceTable <- abundanceTable[,match(samples$Sample_ID,colnames(abundanceTable))]
  colnames(abundanceTable) <- samples$Sample_ID
  abundanceTable <- cbind(abundanceTable,groupTotal=rowSums(abundanceTable, na.rm=TRUE))
  abundanceTable <- rbind(abundanceTable,ReadTotals=colSums(abundanceTable,na.rm=TRUE))
  abundanceTable <- cbind(Level=abtable$level[match(rownames(abundanceTable),abtable$assign)],abundanceTable)
  abundanceTable <- data.frame(Taxon_Name = rownames(abundanceTable),abundanceTable)
  
  dir.create(path=file.path(dir,project),recursive=T,showWarnings=F)
  write.table(abundanceTable,file=file.path(dir,project,outfile),sep="\t",quote=F,col.names=T,row.names=F)
}

### computes abundance tables, based on project
### uses sffinfo and sffile to extract and output only 
"output.genus.reads" <- function(con,project,genus="Lactobacillus"){

  reads <- dbGetQuery(con,
                      paste("SELECT read_data.lucyUnique, read_data.Acc, read_data.lucyLC, read_data.lucyRC, tmpTable.Sample_ID, read_data.Run, read_data.keep
                             FROM read_data, (SELECT pool_metadata.Sample_ID, pool_metadata.Pool, pool_metadata.Reverse_Primer, pool_mapping.Run FROM pool_metadata, pool_mapping
                             WHERE pool_metadata.Project='",project,"'
                             AND pool_metadata.Pool=pool_mapping.Pool) as tmpTable, 
                             (SELECT rdp_report.lucyUnique FROM rdp_report
                              WHERE rdp_report.genus_name = '", genus, "') as rdpTable
                             WHERE tmpTable.Reverse_Primer=read_data.Primer_Code
                             AND tmpTable.Run=read_data.Run
                             AND rdpTable.lucyUnique = read_data.lucyUnique                              
                             GROUP BY read_data.lucyUnique",sep=""))
  
  reads <- reads[reads$read_data.keep == 1,]

  write.table(reads[,2:4],file=paste("454Reads.",genus,".reads",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  sfffiles <-  paste(file.path("Amplicon_SFFfiles_AmpProc",paste(unique(reads$read_data.Run),".sff",sep="")),collapse=" ")
  system(paste("sfffile -i ",paste("454Reads.",genus,".reads",sep=""),
                      " -tr ",paste("454Reads.",genus,".reads",sep=""),
                      " -o ", paste("454Reads.",genus,".sff",sep="")," ", sfffiles, sep=""))
  system(paste("sffinfo -s ",paste("454Reads.",genus,".sff",sep="")," > ",paste("454Reads.",genus,".fasta",sep="")))
#  file.remove(c(paste("454Reads.",genus,".reads",sep=""),paste("454Reads.",genus,".sff",sep="")))
  cat(paste("Done sequences are available in ",paste("454Reads.",genus,".fasta",sep=""),"\n"))
  
}

