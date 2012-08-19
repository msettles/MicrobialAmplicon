########## GENERATE SQLite DB

library(RSQLite)

drv <- SQLite()
con <- dbConnect(drv, dbname="amplicondata.sqlite")

dbGetQuery(con, "
CREATE TABLE read_data (
  Acc CHAR(14) PRIMARY KEY,             -- Accession ID of the Read
  Run CHAR(9) NOT NULL,
  Sample_ID VARCHAR(80) NOT NULL,
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
  version CHAR(3) NOT NULL
);
CREATE INDEX Iacc ON read_data (Acc);
CREATE INDEX Isample_ID ON read_data (Sample_ID);
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
  genus_bootstrap
);
CREATE INDEX IqueryName ON rdp_report (QueryName);
")

dbDisconnect(con)
