
ExtractTimeseries <- function(Instr, Data){
  # Function which extracts the time series of a specific Instrument for 
  #
  # Args:
  #    Instr: A character which needs to be exactly one of the Instruments factor levels. 
  #    Data: A data frame from which to create a time series from. Needs to be on long format with the following column variables
  #          - lane_num
  #          - read_num
  #          - Instrument
  #
  # Returns:
  #    A list where each element of the list is a data.frame object corresponding to measurements of a specific read and lane
  
  read.num.series <- list()
  
  numb_Lanes <- unique(Data$lane_num) %>% length()
  numb_Reads <- unique(Data$read_num) %>% length()
  
  for(Lane in 1:numb_Lanes) {
    for(Read in 1:numb_Reads) {
      df <- filter(Data, lane_num==Lane, read_num==Read, Instrument==Instr)  
      # Unite lane_num and read_num to a 
      df <- unite_(df,"lane_read", c("lane_num","read_num"))
      read.num.series <- c(read.num.series,list(df))
    }
  }
  
  return(read.num.series)
}
 #All.df.reduced
 #First_test <- ExtractTimeseries("HiSeq 3", All.df.reduced)

extractVars <- function(what_variables,list_containing_dataframes){
  # Function which extracts the variables of interest from a list where each element of the 
  # list is a data frame. It is assumed that the variable names in what_variables are contained in each data.frame of the list data
  #
  # Args:
  #     what_variables: A vector with variable names. These variable names needs in 
  #                     each and every data frame of argument "list_containing_dataframes"
  #     list_containing_dataframes: A list where each element is a data frame with atleast column "what_variables".
  #     
  # Returns:
  #     A list where each element is a data frame containing the variables which are listed in argument "what_variables"
  
  lapply(list_containing_dataframes, function(y) y[,what_variables])
}
#Second_test <- extractVars(c("Date","flowcell_id","mean_q"),First_test)

JoinTablesList <- function(list_containing_dataframes,what_to_join_by){
  # Function which joins the elements of a list where each element are data frames containing a variable what_to_join_by
  #
  # Args:
  #     list_containing_dataframes: A list which elements are data frames and where each element contains the variable "what_to_join_by"
  #     what_to_join_by: A variable name which needs to be in each element of the list "list_containing_dataframes"
  
  
  n <- length(list_containing_dataframes)
  joined.table <- left_join(list_containing_dataframes[[1]],list_containing_dataframes[[2]], by=what_to_join_by)
  # continue to join data frames
  for(i in 3:n){
    joined.table <- left_join(joined.table,list_containing_dataframes[[i]],by=what_to_join_by)
  }
  return(joined.table)
}

ExtractBotLvlTimeseries <- function(variable, machine, Data){
  # Function which takes data on long format and places it on wide format, a time series format. 
  
  # Args:
  #     variable: What variables would you like to extract to long format. Needs to be one of the variable names in the Data argument.
  #     machine: What machine would you like to extract data for. Needs to be one of the levels in a Instrument variable
  #     Data: A data frame which should be on long format where the following variables needs to be in the data frame
  #           flowcell_id, cycles, Date, lane_num, read_num
  #           Example: of the Data argument 
  #  
  #       flowcell_id lane_num read_num raw_density pf_density error_rate cycles pct_q30 mean_q     Date      Instrument
  #           (chr)   (fctr)    (int)       (dbl)      (dbl)      (dbl)  (int)   (dbl)    (dbl)    (date)       (fctr)    
  # 1    C13HJACXX        1        1    749552.1   708891.9       0.34    100 0.9478176 37.01983  2012-09-04    HiSeq 3       
  # 2    C13HJACXX        1        2    749552.1   708891.9       0.41    100 0.9187798 36.16825  2012-09-04    HiSeq 3        
  # 3    C13HJACXX        2        1    679662.1   646595.8       0.30    100 0.9540000 37.21812  2012-09-04    HiSeq 3        
  # 8    C13HJACXX        4        2    687972.9   654434.4       0.34    100 0.9293122 36.47404  2012-09-04    HiSeq 3        
  # 9    C13HJACXX        5        1    317835.5   311095.9       0.25    100 0.9798837 38.15038  2012-09-04    HiSeq 3        
  # 10   C13HJACXX        5        2    317835.5   311095.9       0.31    100 0.9707265 37.86044  2012-09-04    HiSeq 3        
  # ........
  #  
  #
  # Returns:
  #         A data frame were each row is a observation and each column corresponds to a observation
  #         for specific variables in a specific lane and read. 
  #         
  
  # join by date and flowcell_id
  by <- c("Date", "flowcell_id", "cycles")
  # what to extract from data
  what <- c("Date","flowcell_id","lane_read","cycles", variable)
  
  # extract all timeseries to a lsit
  tmp.df <- ExtractTimeseries(Instr=machine, Data=Data)
  # Extract the timeseries that we are interested in.
  tmp.df <- extractVars(what,tmp.df)
  # fix name to correspond to which lane and read num
  tmp.df <- lapply(tmp.df, 
                   function(x){
                     col.name <- colnames(x)
                     col <- colnames(x) %in% variable
                     colnames(x)[col] <- paste0(variable,"_",unique(x$lane_read))
                     x$lane_read <- NULL
                     return(x)
                   }
                  )
  
  tmp.df <- JoinTablesList(tmp.df,by) 
  
  return(tmp.df)
}
