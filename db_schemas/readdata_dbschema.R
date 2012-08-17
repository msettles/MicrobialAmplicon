########## GENERATE SQLite DB

library(RSQLite)

drv <- SQLite()
con <- dbConnect(drv, dbname="readdata.sqlite")

dbGetQuery(con, "
CREATE TABLE readdata (
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
  Barcode VARCHAR(80),
  FPErr INTEGER,
  Code_Dist INTEGER,
  Primer_Reverse VARCHAR(80),
  RPErr INTEGER,
  keepAdapter VARCHAR(5) NOT NULL,
  lucyLC INTEGER,
  lucyRC INTEGER,
  lucyLength INTEGER,
  LucyUnique VARCHAR(80) NOT NULL,
  lucyNs INTEGER,
  lucymHomoPrun INTEGER,
  LucyFlip VARCHAR(5) NOT NULL,
  LucyRDPgenus VARCHAR(80),
  LucyRDPboot NUMERIC,
  LucyTE INTEGER,
  LucyQM INTEGER,
  keepLucy VARCHAR(5) NOT NULL
);
CREATE INDEX Iacc ON readdata (Acc);
CREATE INDEX Isample_ID ON readdata (Sample_ID);
CREATE INDEX Irun ON readdata (Run);
CREATE INDEX Ibarcode ON readdata (Barcode);
")

con <- dbConnect(drv, dbname="readdata.sqlite")

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
  SimBtwnQuery.Template NUMERIC,
  QueryFull INTEGER,
  flip VARCHER(5) NOT NULL,
  adpLC INTEGER
};
CREATE INDEX IqueryName ON align_report (QueryName);
")




dbDisconnect(con)
