
Extract.timeseries <- function(x, Data){
  # x: (chr) what machine do you want to use?
  # Data: (data.frame) data which to extract from.
  read.num.series <- list(NULL)
  
  if ( any(colnames(Data) %in% c("read_num")) ){
    n1 <- unique(Data$lane_num) %>% length()
    n2 <- unique(Data$read_num) %>% length()
    
    for(i in 1:n1) {
      for(j in 1:n2) {
        df <- filter(Data, lane_num==i, read_num==j, Instrument==x)  
        df <- unite_(df,"lane_read", c("lane_num","read_num"))
        read.num.series <- c(read.num.series,list(df))
      }
    }
  }else{
    n1 <- unique(Data$lane_num) %>% length()
    
    for(i in 1:n1) {
      df <- filter(Data, lane_num==i, Instrument==x)  
      read.num.series <- c(read.num.series,list(df))
    }
  }
  
  read.num.series[[1]] <- NULL
  return(read.num.series)
}


extractVars <- function(x,data){
  # data: (data.frame) a list of data frames.
  # x: (vector, chr) what to extract from data 
  lapply(data, function(y) y[,x])
}
#extractVars(c("Date","flowcell_id","mean_q"),Extract.timeseries("HiSeq 3", All.df.reduced))

join.tables.list <- function(x,y){
  n <- length(x)
  joined.table <- left_join(x[[1]],x[[2]], by=y)
  
  for(i in 3:n){
    joined.table <- left_join(joined.table,x[[i]],by=y)
  }
  return(joined.table)
}

ExtractBotLvlTimeseries <- function(variable, machine, Data){
  
  # Args:
  # variable: (Vector chr) vector with variable names
  # Data: (data.frame) a dataframe containing the variables lane_num,read_num 
  # machine: (chr) which machine
  
  # Returns:
  # A df in long format consisting of different timeseries in each column

  #####
  
  # join by date and flowcell_id
  by <- c("Date", "flowcell_id", "cycles")
  # what to extract from data
  if (any(colnames(Data) %in% c("read_num")) ){
    what <- c("Date","flowcell_id","lane_read","cycles", variable)
  }else{
    what <- c("Date","flowcell_id","lane_num","cycles", variable)
  }
  # extract all timeseries to a lsit
  tmp.df <- Extract.timeseries(machine, Data)
  # Extract the timeseries that we are interested in.
  tmp.df <- extractVars(what,tmp.df)
  # fix name to correspond to which lane and read num
  if (any(colnames(Data) %in% c("read_num")) ){
    tmp.df <- lapply(tmp.df, function(x){
      col.name <- colnames(x)
      col <- colnames(x) %in% variable
      colnames(x)[col] <- paste0(variable,"_",unique(x$lane_read))
      x$lane_read <- NULL
      return(x)
    })
  }else{
    tmp.df <- lapply(tmp.df, function(x){
      col.name <- colnames(x)
      col <- colnames(x) %in% variable
      colnames(x)[col] <- paste0(variable,"_",unique(x$lane_num))
      x$lane_num <- NULL
      return(x)
    })
  }
  
  
  tmp.df <- join.tables.list(tmp.df,by) 
  
  return(tmp.df)
}
