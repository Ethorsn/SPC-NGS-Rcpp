library(RMySQL)
library(dplyr)
library(lubridate)

# Connect to LOCAL MySql database
PASSWORD <- "A password!"
con <- dbConnect(RMySQL::MySQL(),dbname="ProjectMan", username="root", password=PASSWORD)

# list tables.
# dbListTables(con)

# select everything from qpcr
dbGetQuery(con,"select * from qpcr2012") %>% head()

getDates.id <- "select runfolder_name,flowcell_id from flowcell_runfolder"

df.Dates.id  <- dbGetQuery(con,getDates.id) %>% tbl_df()

# Function which parse out dates from runfolder name.
# Note! This variable is originally contained in the Projman database but was lost in the transition to MySQL.
parseOutDate <- function(x){
  # Parse out the date from the runfolder name variable.
  #
  # Args:
  #   x: A data frame containing a column called "runfolder_name"
  #       on the standard Illumina format: "140314_D00458_0004_BH8CPTADXX",
  #       where the first element is the date.
  #
  # Returns:
  #   The same data frame with a "Date" column with the date on format %y%m%d.
  tmp.split <- lapply(strsplit(x$runfolder_name,"_"), function(x) x[1]) 
  
  tmp.split <- as.Date(unlist(tmp.split), "%y%m%d")
  
  x$Date <- tmp.split
  return(x)
}

###### From data.R Rdashboard.
parseOutInstrumentName <- function(df) {    
  # Parses the instrument name from the runfolder name and
  # translates it to the "common names" of the instruments,
  # e.g. SN334 = HiSeq1
  #
  # Args:
  #   df: A data frame containing a column called "runfolder_name"
  #       on the standard Illumina format: "140314_D00458_0004_BH8CPTADXX",
  #       where the second element is the instrument identifer.
  #
  # Returns:
  #   The same data frame with a "Instrument" column with the instrument
  #   common names.
  
  instrument.translation.table <-
    data.frame(UniqueIdentifier=c("SN344", "SN866", "SN7001335", "D00118", "D00457", "D00458", "ST-E00215", "ST-E00216", "ST-E00274", "ST-E00279", "ST-E00280","M00485", "M00629"),
               Instrument=c("HiSeq 1","HiSeq 2", "HiSeq 3", "HiSeq 4", "HiSeq 5", "HiSeq 6", "HiSeqX 1", "HiSeqX 2", "HiSeqX 3", "HiSeqX 4", "HiSeqX 5","MiSeq 1", "MiSeq 2"))
  
  ExtractSecondElement <- function(x) lapply(x, function(x) x[2])
  
  tmp.list <- strsplit(as.character(df.Dates.id$runfolder_name),"_")
  
  tmp.uniq.identifiers <- data.frame(UniqueIdentifier = do.call(rbind, ExtractSecondElement(tmp.list)))
  
  instrument <- Map(x=tmp.uniq.identifiers$UniqueIdentifier, function(x) {
    instrument.translation.table[which(instrument.translation.table$UniqueIdentifier %in% x),]$Instrument
  })
  
  df$Instrument <- unlist(instrument)
  
  df
}

df.Dates.id <- parseOutDate(df.Dates.id)

df.Dates.id <- parseOutInstrumentName(df.Dates.id) 


Query_flowcell <- "select * from flowcell_lane_results"
df.flowcell <- dbGetQuery(con, Query_flowcell)
df.flowcell <- df.flowcell  %>%  
  mutate(lane_num = as.factor(lane_num),
         read_num = as.integer(read_num)) %>%
  tbl_df()

df.flowcell$delivered <- NULL

######################################## 
# Join tables.
All.df <- left_join(df.flowcell,df.Dates.id, by="flowcell_id") 

All.df <- All.df %>%
  group_by(Instrument) %>% 
  arrange(Date) %>% 
  mutate(DaysIdle=Date-lag(Date))

rm(df.flowcell)
# Remove instruments which are not in use 
All.df.reduced <- All.df[!All.df$Instrument %in% c("HiSeq 1","HiSeq 2","MiSeq 2"),]

# All are NA on these rows.
All.df.reduced[which(is.na(All.df.reduced)),]
# How many rows?
All.nas <- nrow(All.df.reduced[which(is.na(All.df.reduced)),]) 

All.df.reduced <- mutate(All.df.reduced, DaysIdle=as.integer(DaysIdle))

#################################################################
# Using the sample_results table 
Query_sample <- "select * from sample_results" 

df.sample.results <- dbGetQuery(con, Query_sample) %>% tbl_df()
# missing data -> all are NA.
# View(df.sample.results[which(is.na(df.sample.results)),])
# Check if ALL are NA. 
# any(is.na(df.sample.results[which(is.na(df.sample.results)),]) != TRUE)
# Remove these
df.sample.results <- na.omit(df.sample.results)

df.sample.results <- df.sample.results %>%
  mutate(sample_name = as.factor(sample_name), 
         tag_seq = as.factor(tag_seq)) 

df.sample.results <- left_join(df.sample.results,df.Dates.id,by="flowcell_id")
df.sample.results <- df.sample.results[!df.sample.results$Instrument %in% c("HiSeq 1","HiSeq 2","MiSeq 2"),]

###################################################################
# disconnect to database
dbDisconnect(con)
# close
rm(df.Dates.id,addClassfun,parseOutDate,parseOutInstrumentName,Query_flowcell,Query_sample,con,getDates.id,All.df)
