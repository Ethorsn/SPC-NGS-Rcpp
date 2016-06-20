## ----RlibsFuns, echo=FALSE, chache=TRUE----------------------------------
# Libraries
library(ggplot2)
library(reshape2)
library(GGally)
library(lubridate)
library(tidyr)
library(dplyr)
library(portes)
library(car) 
library(gridExtra)
library(grid)
library(ggExtra)
library(parallel)
library(xtable)
# load data + functions
# must be run in this order!!!
source("../FunctionsAndRcpp/VizCorFuns.R")
#source("../FunctionsAndRcpp/mysqlScript.R")
load("../Data/MySQL_data.Rdata")
source("../FunctionsAndRcpp/BotlvlSeriesExtraction.R")

###################### Hadley Wickhams function
grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}
#######################

## ----Taglvl--------------------------------------------------------------
N <- unique(df.sample.results$flowcell_id) %>% length()
First.date <- first(df.sample.results$Date)
last.date <- last(df.sample.results$Date)

## ----FilterSample--------------------------------------------------------
# Number of unique flowcells
N2 <- unique(All.df.reduced$flowcell_id) %>% length()
Missing <- anti_join(All.df.reduced,df.sample.results, by="flowcell_id") %>% .$flowcell_id %>% unique() %>% length()

# exclude data before 2013.
df.sample.results <- filter(df.sample.results, year(Date)>2012)
All.df.reduced <- filter(All.df.reduced,year(Date)>2012)


## ----Table, results='asis', warnings=FALSE, eval=FALSE-------------------
## # These transformations are to create the table containing the number of flowcells on a specific completed run cycle.
## ggplot.tmp <- df.sample.results[df.sample.results$Instrument != "MiSeq 1",] %>%
##   group_by(Instrument,cycles) %>%
##   summarise(x=length(unique(flowcell_id)))
## 
## ggplot.tmp$cycles <- as.factor(ggplot.tmp$cycles)
## 
## xtab <- spread(ggplot.tmp,key="Instrument",value=x)
## colnames(xtab) <- c("cycles", paste0("Hi",3:6), paste0("HiX",1:5))
## xtab[is.na(xtab)] <- 0
## 
## xtab <- rbind(xtab,c(NA,apply(xtab[,-1],2,sum,na.rm=TRUE)))
## 
## print(xtable(xtab,
##        caption = "Table showing the number of completed cycles for each HiSeq (Hi) and HiSeqX (HiX) machine. Notice that the MiSeq 1 machine is removed.",
##        label="CompCycl",
##        digits=0),
##   include.rownames=FALSE)
## # The table presented in the Thesis is somewhat modified in the appearance but NOT in the content.

## ----TagLevelTS, fig.cap="Figure containing the range (min to max) and mean of each successive run (flowcell). Here, we are showing read 1 and 2 in lane 1, disregarding what type of setting the run is performed on.", fig.height=7,fig.width=10, fig.pos="!htbp"----

# Create data.frame containing mean and range of mean q measurements on tag level
ggplot.tmp1 <- df.sample.results %>%
  filter(Instrument == "HiSeq 6" |
         Instrument == "MiSeq 1" |
         Instrument == "HiSeqX 1", 
         lane_num == 1) %>%
  group_by(Date,flowcell_id, Instrument, read_num) %>%
  summarise(mn = mean(mean_q), min = min(mean_q), max = max(mean_q)) %>%
  group_by() %>% 
  arrange(desc(Instrument)) 
# How many observations for each machine?
n1 <- filter(ggplot.tmp1, Instrument == "HiSeq 6") %>% nrow()
n2 <- filter(ggplot.tmp1, Instrument == "HiSeqX 1") %>% nrow()
n3 <- filter(ggplot.tmp1, Instrument == "MiSeq 1") %>% nrow()

# Order of appearance: MiSeq 1, HiSeqX 1, HiSeq 6  
OrderOfApp <- c(rep(seq(1,n3/2, by=1), each=2), rep(seq(1,n2/2, by=1),each=2), rep(seq(1,n1/2, by=1),each=2))
# Mutate Date to order of appearance.
ggplot.tmp1 <- mutate(ggplot.tmp1, Date = OrderOfApp)
# Rename variables for nicer printing
ggplot.tmp1 <- rename(ggplot.tmp1, Machine=Instrument, Read=read_num)
# Create figure of mean and range of mean q measurements on tag level
ggplot(ggplot.tmp1,
       aes(x=Date,y=mn)) +
  geom_errorbar(aes(ymin=min,
                    ymax=max,
                    color="#7fcdbb"),
                stat="identity",
                width=0.5,
                alpha=0.8) +
  scale_shape_identity() +
  geom_point(aes(color="#990000"), alpha=0.6) +
  facet_grid(Read~Machine,
             scales = "free",
             #ncol=2,
             labeller = label_both) +
  theme_bw() +
  ylab("value") +
  scale_color_manual("",
                     label=c("Range (min-max)","mean"),
                     values=c("#7fcdbb","#990000")) +
  theme(axis.text.x=element_text(angle=60,
                                 vjust=0.7,
                                 hjust=0.8),
        legend.position="bottom") +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ggtitle("mean q score measurements on Tag level") +
  xlab("Order of appearance")

## ----HiSeq6Comb, fig.cap="Figure showing the range (min to max) and mean of each succsessive run (flowcell) in lane 1, read 1, of the two variables Percent q30 and Percentage tag error. All runs shown where performed on 126 cycles.", fig.pos="!htb"----

# Completed run cycle. Assumed to represent the actual run cycles
Cycl <- c(102,126)
# create a data.frame for plotting the range and mean of variables percent q30 and percent tag error, lane 1 read 1. 
ggplot.tmp <- df.sample.results %>%
  filter(Instrument == "HiSeq 6",
         cycles>=Cycl[1],
         cycles<=Cycl[2],
         lane_num == 1,
         read_num == 1) %>%
  group_by(flowcell_id, Date) %>%
  summarise(mnPQ30 = mean(pct_q30/100), minPQ30 = min(pct_q30/100), maxPQ30 = max(pct_q30/100),
            mnPTerr = mean(pct_tag_err/100), minPTerr = min(pct_tag_err/100), maxPTerr = max(pct_tag_err/100),
            NumbObs=length(pct_q30/100)) %>%
  arrange(desc(Date))

# Figure containing mean and range of percent q30 variable for read = 1, lane = 1
p1 <- ggplot(ggplot.tmp,
       aes(x=seq(0,length(Date)-1),y=mnPQ30)) +
  geom_errorbar(aes(ymin=minPQ30,
                    ymax=maxPQ30,
                    color="#7fcdbb"),
                stat="identity",
                width=0.5,
                alpha=0.8) +
  scale_shape_identity() +
  geom_point(aes(color="#990000"), alpha=0.6) +
  theme_bw() +
  ylab("value") +
  scale_color_manual("",
                     label=c("Range (min-max)","mean"),
                     values=c("#7fcdbb","#990000")) +
  ggtitle("Percent q30") +
  theme(axis.text.x=element_text(angle=60,
                                 vjust=0.7,
                                 hjust=0.8),
        axis.title.y=element_text(size=8),
        title = element_text(size=9),
        legend.position="bottom") +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  xlab("") +
  scale_y_continuous(limits=c(0,1)) 

# Figure containing mean and range of percent tag error variable for read = 1, lane = 1
p2 <- ggplot(ggplot.tmp,
       aes(x=seq(0,length(Date)-1),y=mnPTerr)) +
  geom_errorbar(aes(ymin=minPTerr,
                    ymax=maxPTerr,
                    color="#7fcdbb"),
                stat="identity",
                width=0.5,
                alpha=0.8) +
  scale_shape_identity() +
  geom_point(aes(color="#990000"), alpha=0.6) +
  theme_bw() +
  ggtitle("Percent tag error") +
  ylab("value") +
  scale_color_manual("",
                     label=c("Range (min-max)","mean"),
                     values=c("#7fcdbb","#990000")) +
  theme(axis.text.x=element_text(angle=60,
                                 vjust=0.7,
                                 hjust=0.8),
        axis.title.y=element_text(size=8),
        title = element_text(size=9),
        legend.position="bottom") +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0)))) +
  xlab("")
# Figure with number of observations in each read = 1, lane = 1. 
p3 <- ggplot(ggplot.tmp,aes(x=seq(0,length(Date)-1),y=NumbObs)) +
  geom_point(alpha=0.8) +
  theme_bw() +
  ylab("value") +
  ggtitle("Number of observations") +
  xlab("Order of appearance") +
  theme(axis.text.x=element_text(angle=60,
                                 vjust=0.7,
                                 hjust=0.8),
        title = element_text(size=9),
        axis.title.y=element_text(size=8)) 
grid_arrange_shared_legend(p1,p2,p3)
rm(p1,p2,p3)

## ----MeanVectorFigure, fig.cap="Mean together with the range of the error rate of each lane and read (lane\\_read). Notice that the HiSeq 5 has the lowest mean error rate of all HiSeq machines.", fig.pos="!htb", fig.height=4----

# Extract HiSeq 4 data on wide format. Remove rapid runs and then select error rate.
Hiseq4 <- ExtractBotLvlTimeseries(c("error_rate"),"HiSeq 4",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  # select error rate variable for each lane and read
  select(grep("error_rate", colnames(.)))
# Create Mean, max, min for each lane and read.
Hiseq4 <- data.frame("Mean"=colMeans(Hiseq4),
                     "Max"=apply(Hiseq4,2, max),
                     "Min"=apply(Hiseq4,2, min),
                     "lane_read"=as.factor(colnames(Hiseq4)),
                     "Machine"=as.factor("HiSeq 4")) %>%
  # place into long format for plotting
  melt(measure.vars=c("Mean"))

# Do the same for HiSeq 5 data
Hiseq5 <- ExtractBotLvlTimeseries(c("error_rate"),"HiSeq 5",All.df.reduced) %>%
    filter(cycles>=Cycl[1],
           cycles<=Cycl[2]) %>%
   na.omit() %>%
   select(grep("error_rate", colnames(.)))

Hiseq5 <- data.frame("Mean"=colMeans(Hiseq5),
                     "Max"=apply(Hiseq5,2, max),
                     "Min"=apply(Hiseq5,2, min),
                     "lane_read"=as.factor(colnames(Hiseq5)),
                     "Machine"=as.factor("HiSeq 5")) %>%
  melt(measure.vars=c("Mean"))

# Do the same for HiSeq 6 data
Hiseq6 <- ExtractBotLvlTimeseries(c("error_rate"),"HiSeq 6",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  select(grep("error_rate", colnames(.)))

Hiseq6 <- data.frame("Mean"=colMeans(Hiseq6),
                     "Max"=apply(Hiseq6,2, max),
                     "Min"=apply(Hiseq6,2, min),
                     "lane_read"=as.factor(colnames(Hiseq6)),
                     "Machine"=as.factor("HiSeq 6")) %>%
  melt(measure.vars=c("Mean"))

ggplot.tmp <- rbind(Hiseq4,Hiseq5,Hiseq6)
levels(ggplot.tmp$lane_read) <- do.call("rbind",strsplit(levels(ggplot.tmp$lane_read),"error_rate_"))[,2] 

# Create mean and range figure for error rate.
ggplot(ggplot.tmp,
       aes(y=value, x=lane_read, color=Machine)) +
  geom_point(aes(fill=variable), size=2, position=position_dodge(width=0.5), alpha=0.8) +
  geom_errorbar(aes(ymax=Max, ymin=Min, width=0.2, fill=as.factor(2)), alpha=0.5, position=position_dodge(width=0.5), show.legend = TRUE) +
  scale_y_continuous(limits=c(-0.01,2), breaks = c(0,0.5,1,1.5,2)) +
  coord_flip() +
  theme_bw() +
  scale_color_brewer(type = "qual", palette=6,guide = guide_legend(reverse=TRUE)) +
  scale_fill_manual("",
                    labels=c("Mean","Range (min-max)"),
                    values =c("black","blue")) +
  guides(fill=guide_legend(override.aes=list(shape=c(16,NA), linetype=c(0,1)))) +
  xlab("Error rate (lane_read)") 

## ----Readlvlseries-------------------------------------------------------
# Extract HiSeq 6 data, all variables on read level
HiSeq6.botlvl.series <- ExtractBotLvlTimeseries(c("mean_q","pct_q30","error_rate","raw_clusters","raw_density","pf_clusters","pf_density"),"HiSeq 6",All.df.reduced)
## Extract "High output v4" runs: 102-126 cycles using 8 lanes
v4.HiSeq6.botlvl.series <- filter(HiSeq6.botlvl.series, 
                                  cycles>=Cycl[1],
                                  cycles<=Cycl[2]) %>% 
  # now remove rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit()

## ----ReadlvlER, fig.cap="Error rates, the raw density and the number of raw clusters for each read in lanes 1 and 2. The variable name are listed in the following manner: Variable\\_lane\\_read.", fig.height=4.5, fig.pos="!htb"----

# Code found online and has been slightly modified.
# http://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_printing <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     for (i in 1:length(l)){
       if (grepl("NA",l[i])==TRUE){
         l[i] <- l[i]
       }else{
         if (as.numeric(l[i])>2){
           # quote the part before the exponent to keep all the digits
           l[i] <- gsub("^(.*)e", "'\\1'e", l[i])
           # turn the 'e+' into plotmath format
           l[i] <- gsub("e", "%.%10^", l[i])
         }else{
          l[i] <- as.numeric(l[i])
         }
       }
     }
     # return this as an expression
     parse(text=l)
}
#fancy_printing(c("   NA","1e8"))
ggplot(melt(
  select(v4.HiSeq6.botlvl.series,
    starts_with("error_rate_1"),
    starts_with("error_rate_2"),
    starts_with("raw_density_1"),
    starts_with("raw_density_2"),
    starts_with("raw_clusters_1"),
    starts_with("raw_clusters_2"),
    Date),
  id.vars="Date"),
  aes(x=value)) +
  geom_histogram() +
  facet_wrap("variable", scales = "free_x",ncol=4) +
  theme_bw() +
  theme(legend.position="none", plot.margin=unit(c(5.5, 15, 5.5, 5.5),"points"))  +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3), labels=fancy_printing)

## ----ReadlvlCor,fig.cap="Spearman correlation matrix of HiSeq 6 Read level measurements. Two groupings can be seen in the correlation.", warning=FALSE----
ggp <- VizCor(v4.HiSeq6.botlvl.series[,-c(1:3)], method="spearman",title="")
N <- ncol(v4.HiSeq6.botlvl.series)

tmpFUN <- function(x) return (cbind(c(2+8*(1+2*x)),c(N-8*(1+2*x))))
df.tmp <- tmpFUN(seq(0,6))

df.tmp <- data.frame(TEXT=c("Error rate","Mean Q","Percent q30","PF clusters","PF density","Raw clusters","Raw density"),df.tmp)

ggp <- ggp + 
   theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) 

# add text to picture
for (i in 1:nrow(df.tmp)){
  ggp <- ggp + 
   geom_text(x=df.tmp[i,2], y=df.tmp[i,3],label=df.tmp[i,1], 
             angle=-45,size=3.2,family="Times") 
}
tmpFUN <- function(x) return(cbind(N-(2+16*x)-1/2,16*x+1/2))
df.tmp <- tmpFUN(1:7) %>% as.data.frame()  
# Draw lines to create sections
for (i in 1:nrow(df.tmp)){
  ggp <- ggp + 
    geom_segment(y=df.tmp[i,1], x=0.5, xend=df.tmp[i,2], yend=df.tmp[i,1], color="#0000004D", size=0.2) +
    geom_segment(x=df.tmp[i,1], y=0.5, xend=df.tmp[i,1], yend=df.tmp[i,2], color="#0000004D", size=0.2) 
}
# print
ggp
# remove
rm(df.tmp, tmpFUN)

## ----ICsampleCreation----------------------------------------------------
# Create a IC sample by removing values runs with large error rates and low pf_clusters
Hiseq6.tag <- 
  filter(df.sample.results,
         Instrument == "HiSeq 6", 
         cycles < Cycl[2],
         cycles > Cycl[1]) 

Hiseq6.read <- filter(All.df.reduced,
                      Instrument=="HiSeq 6", 
                      cycles < Cycl[2],
                      cycles > Cycl[1])

CreateICsample <- function(x, y,quality_crit){
  # Function which creates a in control sample for tag and read level. 
  #
  # Args:
  #     x: A data frame which contains "tag level" measurements. 
  #        It is assumed that the variables specified in quality_crit are contained in this data frame.
  #     y: A data frame which contains "read level" measurements.
  #        It is assumed that the variables specified in quality_crit are contained in this data frame.
  #     quality_crit: A list containing the variables and limits which we would like to use as quality limits.
  #                   A example would be quality_crit=list(error_rate=1, pct_q30=0.7).
  #
  # Returns:
  #     Returns a list where each element is a data frame representing a in control data frame. 
  #     The list will have two element corresponding to tag and read level.
  library(dplyr)
  IC.read <- y %>%
    filter(error_rate < quality_crit$err,
           pf_clusters > quality_crit$pf_clusters) 
  
  # How many flowcells are removed?
  # FlowcellRemoved <- anti_join(IC.read,y,by="flowcell_id")$flowcell_id %>% unique() 
  
  # Remove flowcells which are not ok to the quality characteristics
  Data <- semi_join(x,IC.read,by="flowcell_id")
  
  return(list(IC.tag = Data,
              IC.read = IC.read))
}

IC.Hiseq6.sample <- CreateICsample(x=Hiseq6.tag, 
                                   y=Hiseq6.read,
                                   quality_crit = list(err=1, 
                                                       pf_clusters=1e8))
# Extract the data from this in control data set.
IC.Hiseq6 <- ExtractBotLvlTimeseries(c("mean_q","pct_q30","error_rate"),
                                     machine="HiSeq 6", 
                                     Data=group_by(IC.Hiseq6.sample$IC.read)) %>%
  na.omit()

rm(df.sample.results,All.df.reduced,HiSeq6.botlvl.series,Hiseq6.tag,Hiseq6.read,Missing)

## ----TransformationIC----------------------------------------------------
TransformDataFun <- function(x){
  # Function which transforms data using a Box-Cox transformation. 
  # 
  # Args:
  #     x: a vector containing numeric data.
  #
  # Returns: 
  #     A list containing two elements. The first element of the list is the estimated transformation parameter. The second is the transformed data.
  
  # The estimation of the transformation parameter.
  lambda <- forecast::BoxCox.lambda(x)
  
  return(list(Lambda=lambda,
              data=forecast::BoxCox(x,lambda)))
}
##############################
# Transforming a in control sample containing the following variables:
# pct_q30, error_rate, mean_q

# pct_q30 will be transformed using normal quantile function
colnam <- colnames(IC.Hiseq6[,-c(1:3)])
PCT <- grep("pct",colnam)

# Prepare data frame.
TransformedData <- IC.Hiseq6[,-c(1:3)]
# Transform using Box-Cox
TransformedData.tmp <- apply(TransformedData[,!1:ncol(TransformedData) %in% PCT], 2, TransformDataFun)

Lambdas <- lapply(TransformedData.tmp,function(x) x$Lambda) %>% unlist()
# Extract Box-Cox transformaed data.
TransformedData.tmp <- lapply(TransformedData.tmp, function(x) x$data) %>% do.call("cbind",.)
# replace non-transformed data with transformed.
TransformedData[,!1:ncol(TransformedData) %in% PCT] <- TransformedData.tmp
# Transform pct_q30 with quantile function
TransformedData[,PCT] <- apply(TransformedData[,1:ncol(TransformedData) %in% PCT],2,function(x) qnorm(x))
#save(Lambdas, TransformedData, file="../Data/ICdata.Rdata")

## ----Tests---------------------------------------------------------------
# Perform statistical tests on the transformed data. Is it normally distributed?
# Place on the same scale -> divide by a constant.
K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

t2 <- MVN::hzTest(TransformedData)
t4 <- mvShapiroTest::mvShapiro.Test(as.matrix(TransformedData))
#save(t2,t4,file = "../Data/LambdaAndXtable.Rdata")

