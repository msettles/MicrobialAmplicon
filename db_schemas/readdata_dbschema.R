########## GENERATE SQLite DB

#### Version 2.0db
library(RSQLite)

drv <- SQLite()
con <- dbConnect(drv, dbname="amplicondataV2.0.sqlite")

dbBeginTransaction(con)

# Processed files
sql <- "CREATE TABLE processed (
  Run CHAR(9) PRIMARY KEY
  Reads_good INTEGER,
  Reads_fail INTEGER
);"
dbSendQuery(con, sql)
sql <- "CREATE INDEX proc_run ON processed (Run);"
dbSendQuery(con, sql)

# Generate table read_data, stores primary information about each read
sql <- "
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
  Barcode VARCHAR(80),
  FPErr INTEGER,
  Code_Dist INTEGER,
  Primer_3prime VARCHAR(80),
  RPErr INTEGER,
  lucyLC INTEGER,
  lucyRC INTEGER,
  lucyLength INTEGER,
  lucyUnique VARCHAR(80) NOT NULL,
  lucyNs INTEGER,
  lucymHomoPrun INTEGER,
  keep VARCHAR(5) NOT NULL,
  version VARCHAR(80) NOT NULL
);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rd_acc ON read_data (Acc);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rd_run ON read_data (Run);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rd_barcode ON read_data (barcode);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rd_lucy_unique ON read_data (lucyUnique);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rd_keep ON read_data (keep);"
dbSendQuery(con, sql)

# Generate table for alignment data, stores information about mothur alignments
sql <- "
CREATE TABLE align_report (
  lucyUnique VARCHAR(80) PRIMARY KEY,
  Run CHAR(9) NOT NULL,
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
);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX ar_lucy_unique ON align_report (lucyUnique);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX ar_run ON align_report (Run);"
dbSendQuery(con, sql)

# Generate table for rdp_report, stores rdp assignment information
sql <- "CREATE TABLE rdp_report (
  lucyUnique VARCHAR(80)  PRIMARY KEY,
  Run CHAR(9) NOT NULL,
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
  species_name
);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rdp_lucy_unique ON rdp_report (lucyUnique);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rdp_run ON rdp_report (Run);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX rdp_genus ON rdp_report (genus_name);" ## index on genus, for speciation
dbSendQuery(con, sql)

# Generate table to store pool metadata
sql <- "CREATE TABLE pool_metadata (
  ID VARCHAR(20) NOT NULL,
  Pool VARCHAR(20) NOT NULL,
  Second_Primer VARCHAR(20) NOT NULL,
  Barcode VARCHAR(80), NOT NULL,
  DNA_conc NUMERIC,
  tVol  NUMERIC,
  pVol NUMERIC,
  Investigator VARCHAR(20) NOT NULL,
  Prepared_by VARCHAR(20) NOT NULL,
  Isolation_ID VARCHAR(20) NOT NULL,
  Project VARCHAR(20) NOT NULL,
  Sample_ID VARCHAR(20) PRIMARY_KEY
);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX pm_pool ON pool_metadata (Pool);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX pm_barcode ON pool_metadata (Barcode);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX pm_project ON pool_metadata (Project);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX pm_sample_id ON pool_metadata (Sample_ID);"
dbSendQuery(con, sql)

# Generate table 
sql <- "
CREATE TABLE pool_mapping (
  Pool VARCHAR(20) NOT NULL,
  Run CHAR(9) NOT NULL
);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX map_pool ON pool_mapping (Pool);"
dbSendQuery(con, sql)

sql <- "CREATE INDEX map_run ON pool_mapping (Run);"
dbSendQuery(con, sql)


dbCommit(con)
dbDisconnect(con)
