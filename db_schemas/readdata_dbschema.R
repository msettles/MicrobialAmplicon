########## GENERATE SQLite DB

library(RSQLite)

drv <- SQLite()
con <- dbConnect(drv, dbname="amplicondata.sqlite")

dbGetQuery(con, "
CREATE TABLE read_data (
  Acc CHAR(14) PRIMARY KEY,             -- Accession ID of the Read
  Run CHAR(9) NOT NULL,
  RawLength INTEGER,
  RocheLC INTEGER,
  RocheRC INTEGER,
  RocheLength INTEGER,
  AdapterLC INTEGER,
  AdapterRC INTEGER,
  AdapterLength INTEGER,
  Primer_Code VARCHAR(80),
  FPErr INTEGER,
  Barcode VARCHAR(80),
  Code_Dist INTEGER,
  Primer_Reverse VARCHAR(80),
  RPErr INTEGER,
  lucyLC INTEGER,
  lucyRC INTEGER,
  lucyLength INTEGER,
  lucyUnique VARCHAR(80) NOT NULL,
  lucyNs INTEGER,
  lucymHomoPrun INTEGER,
  keep VARCHAR(5) NOT NULL,
  version VARCHAR(80) NOT NULL
);
CREATE INDEX Iacc ON read_data (Acc);
CREATE INDEX Irun ON read_data (Run);
CREATE INDEX Ibarcode ON read_data (Barcode);
")

dbGetQuery(con, "
CREATE TABLE align_report (
  QueryName VARCHAR(80) PRIMARY KEY,
  QueryLength INTEGER,
  TemplateName VARCHAR(20),
  TemplateLength INTEGER,
  SearchMethod CHAR(4) NOT NULL,
  SearchScore INTEGER,
  AlignmentMethod CHAR(9) NOT NULL,
  QueryStart INTEGER,
  QueryEnd INTEGER,
  TemplateStart INTEGER,
  TemplateEnd INTEGER,
  PairwiseAlignmentLength INTEGER,
  GapsInQuery INTEGER,
  GapsInTemplate INTEGER,
  LongestInsert INTEGER,
  SimBtwnQuery_Template NUMERIC,
  flip VARCHER(5) NOT NULL,
  QueryFull INTEGER,
  adpLC INTEGER
);
CREATE INDEX IqueryName ON align_report (QueryName);
")


dbGetQuery(con, "
CREATE TABLE rdp_report (
  QueryName VARCHAR(80) PRIMARY KEY,
  flip VARCHER(5),
  domain_name,
  domain_bootstrap,
  phylum_name,
  phylum_bootstrap,
  class_name,
  class_boostrap,
  order_name,
  order_boostrap,
  family_name,
  family_bootstrap,
  genus_name,
  genus_bootstrap,
  species_name,
  species_bootstrap
);
CREATE INDEX IqueryName ON rdp_report (QueryName);
")


dbGetQuery(con,"
CREATE TABLE pool_metadata (
  ID VARCHAR(20) NOT NULL,
  Pool VARCHAR(20) NOT NULL,
  Forward_Primer VARCHAR(20) NOT NULL,
  Reverse_Primer VARCHAR(20) NOT NULL,
  DNA_conc NUMERIC,
  tVol  NUMERIC,
  pVol NUMERIC,
  Investigator VARCHAR(20) NOT NULL,
  Prepared_by VARCHAR(20) NOT NULL,
  Isolation_ID VARCHAR(20) NOT NULL,
  Project VARCHAR(20) NOT NULL,
  Sample_ID VARCHAR(20) NOT NULL
);
CREATE INDEX Ppool ON pool_metadata (Pool);
CREATE INDEX Preverse_primer ON pool_metadata (Reverse_Primer);
CREATE INDEX Pproject ON pool_metadata (Project);
CREATE INDEX Psample_id ON pool_metadata (Sample_ID);
")

dbGetQuery(con,"
CREATE TABLE pool_mapping (
  Pool VARCHAR(20) NOT NULL,
  Run CHAR(9) NOT NULL
);
CREATE INDEX Mpool ON pool_mapping (Pool);
CREATE INDEX Mrun ON pool_mapping (Run);
")

poolD <- read.table("MetaData/Pool_MetaData.txt",sep="\t",header=T)
dbWriteTable(con,"pool_metadata",poolD,row.names=F,append=T)

poolM <- read.table("MetaData/Pool_Mapping.txt",sep="\t")
dbWriteTable(con,"pool_mapping",poolM,row.names=F,append=T)

#getsample <-dbGetQuery(con,"Select Run, pool_metadata.Pool, Reverse_Primer, Sample_ID from pool_metadata, pool_mapping WHERE pool_metadata.Pool=pool_mapping.Pool AND pool_metadata.project='Adolescence'")
#getreads <- dbGetQuery(con,"Select pool_metadata.Sample_ID, read_data.* FROM pool_metadata, read_data, pool_mapping WHERE pool_metadata.project='Adolescence' AND pool_metadata.Pool=pool_mapping.Pool AND pool_metadata.Reverse_Primer=read_data.Primer_Code AND pool_mapping.Run=read_data.Run")

dbDisconnect(con)
