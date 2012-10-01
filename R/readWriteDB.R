### Retreive and Send metadata or mapping data


library(RSQLite)

drv <- SQLite()
con <- dbConnect(drv, dbname="amplicondata.sqlite")
## get the current table
mappings <- dbReadTable(con,"pool_mapping")
pooltable <- read.table("MetaData/Pool_Mapping.txt",sep="\t")
pooltable <- pooltable[!(apply(pooltable,1,paste,collapse=" ") %in% apply(mappings,1,paste,collapse=" ")),]
print(paste("Adding ",nrow(pooltable)," entries",sep=""))
if (nrow(pooltable) > 0)
  dbWriteTable(con,"pool_mapping",pooltable,row.names=F,append=T)

#dbGetQuery(con, "DELETE FROM pool_mapping;")
#poolD <- read.table("MetaData/Pool_MetaData.txt",sep="\t",header=T)
#poolD <- read.table("MetaData/Pool_MetaData2.txt",sep="\t",header=T)
#poolD <- read.table("MetaData/Pool_MetaData3.txt",sep="\t",header=T)
#dbWriteTable(con,"pool_metadata",poolD,row.names=F,append=T)


#getsample <-dbGetQuery(con,"Select Run, pool_metadata.Pool, Reverse_Primer, Sample_ID from pool_metadata, pool_mapping WHERE pool_metadata.Pool=pool_mapping.Pool AND pool_metadata.project='Adolescence'")
#getreads <- dbGetQuery(con,"Select pool_metadata.Sample_ID, read_data.* FROM pool_metadata, read_data, pool_mapping WHERE pool_metadata.project='Adolescence' AND pool_metadata.Pool=pool_mapping.Pool AND pool_metadata.Reverse_Primer=read_data.Primer_Code AND pool_mapping.Run=read_data.Run")

library(RSQLite)

## Metadata Table to add
mdTable <- "MetaData/Pool_MetaData_Top.txt"
drv <- SQLite()
con <- dbConnect(drv, dbname="amplicondata.sqlite")
## get the current table
metadata <- dbReadTable(con,"pool_metadata")
metatable <- read.table(mdTable,sep="\t",header=T)
metatable <- metatable[!(apply(metatable[,c("Project","Sample_ID","Reverse_Primer")],1,paste,collapse=" ") %in% apply(metadata[,c("Project","Sample_ID","Reverse_Primer")],1,paste,collapse=" ")),]
print(paste("Adding ",nrow(metatable)," entries",sep=""))
if (nrow(pooltable) > 0)
  dbWriteTable(con,"pool_metadata",metatable,row.names=F,append=T)


getsample <-dbGetQuery(con,"Select Run, pool_metadata.Pool, Reverse_Primer, Sample_ID from pool_metadata, pool_mapping WHERE pool_metadata.Pool=pool_mapping.Pool AND pool_metadata.project='JJ_Human_Vagina'")
getreads <- dbGetQuery(con,"Select pool_metadata.Sample_ID, read_data.* FROM pool_metadata, read_data, pool_mapping WHERE pool_metadata.project='JJ_Human_Vagina' AND pool_metadata.Pool=pool_mapping.Pool AND pool_metadata.Reverse_Primer=read_data.Primer_Code AND pool_mapping.Run=read_data.Run")

