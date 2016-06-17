library(RMySQL)
library(dplyr)
library(lubridate)

con <- dbConnect(RMySQL::MySQL(),dbname="ProjectMan", username="root", password="Mkonji123")
# summary(con)
# list tables.
dbListTables(con)

dbGetQuery(con,"select * from qpcr2012") %>% head()

getDates.id <- "select runfolder_name,flowcell_id from flowcell_runfolder"
df.Dates.id  <- dbGetQuery(con,getDates.id) %>% tbl_df()
#df.Dates.id <- fetch(Dates.id) %>% tbl_df()
#dbClearResult(Dates.id)
parseOutDate <- function(x){
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

parseOutRunType <- function(x) {
  # x: (data.frame) 
  x<- df.Dates.id
  if (!"flowcell_id" %in% colnames(x)){
    print("flowcell id is not contained in this table") 
  }else{
    x$RunType <- ""
    # HiSeq machines Runtypes 
    x$RunType[grep("ANXX", x$flowcell_id)] <- "v4"
    x$RunType[grep("BCXX", x$flowcell_id)] <- "Rapid"
  } 
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
All.df <- All.df[!All.df$Instrument %in% c("HiSeq 1","HiSeq 2","MiSeq 2"),]
# Remove those runs which has 0 cycles.
# All.df <- All.df[!All.df$cycles==0,]
# filter our years before 2013, as Johan suggested.
All.df.reduced <- filter(All.df, year(Date) > 2012)

addClassfun <- function(x){
  Class.instrument <- lapply(strsplit(as.character(x$Instrument)," "), function(x) x[1])  
  Class.instrument<- unlist(Class.instrument)
  x$class <- as.factor(Class.instrument)
  return(x)
}

All.df.reduced <- addClassfun(All.df.reduced)
# All are NA on these rows.
All.df.reduced[which(is.na(All.df.reduced)),]
# How many?
All.nas <- nrow(All.df.reduced[which(is.na(All.df.reduced)),]) 
# Omit.
All.df.reduced <- na.omit(All.df.reduced) # remove them.
# 
All.df.reduced <- mutate(All.df.reduced, DaysIdle=as.integer(DaysIdle))

# df.lane.level <- All.df.reduced %>% 
#   group_by(Instrument, flowcell_id, Date, lane_num,cycles) %>%
#   summarise_each(funs(mean(.,na.rm=TRUE)),
#                  raw_density,error_rate,raw_clusters, raw_clusters, pf_clusters,mean_q)
# 
# df.lane.level <- addClassfun(df.lane.level)
# 
# df.top.level <- df.lane.level %>%
#   group_by(Instrument,flowcell_id,Date) %>%
#   summarise_each(funs(mean(.,na.rm=TRUE)),
#                  raw_density,error_rate,raw_clusters, raw_clusters, pf_clusters,mean_q)
#   
# df.top.level <- addClassfun(df.top.level)

# All.df.scaled[,numerics] <- data.frame(lapply(All.df.reduced[,numerics],
#                                               scale, # function 
#                                               center=TRUE)) # arguments to scale fun.

#################################################################
# Using the sample_results table 
Query_sample <- "select * from sample_results" 
# no date variable -> no ordering.
# OBSOBSOBS this is on sample level, i.e. for each sample.
df.sample.results <- dbGetQuery(con, Query_sample) %>% tbl_df()

# missing data. if mean_q is missing <-> all others are missing
df.sample.results <- na.omit(df.sample.results)

df.sample.results <- df.sample.results %>%
  mutate(sample_name = as.factor(sample_name), 
         tag_seq = as.factor(tag_seq)) 

df.sample.results <- left_join(df.sample.results,df.Dates.id,by="flowcell_id")
df.sample.results <- df.sample.results[!df.sample.results$Instrument %in% c("MiSeq 2","HiSeq 1", "HiSeq 2"),]
##### Do this in MySQL instead. 
rm(df.Dates.id)
dbDisconnect(con)

rm(addClassfun,parseOutDate,parseOutRunType,parseOutInstrumentName,Query_flowcell,Query_sample,con,getDates.id)
