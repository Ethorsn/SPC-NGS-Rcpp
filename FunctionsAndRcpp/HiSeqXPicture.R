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
library(Rcpp)
library(RcppArmadillo)
# load data + functions
# must be run in this order!!!
source("../FunctionsAndRcpp/VizCorFuns.R")
source("../FunctionsAndRcpp/mysqlScript.R")
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
N <- unique(All.df$flowcell_id) %>% length()
Numb2012 <- filter(All.df, year(Date) == 2012) %>% select(flowcell_id) %>% unique() %>% nrow()
First.date <- first(df.sample.results$Date)
###########################################################
Cycl <- c(140,160)
meanQ_HiseqX1 <- ExtractBotLvlTimeseries(c("mean_q"),"HiSeqX 1",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  select(grep("mean", colnames(.)))

meanQ_HiseqX1 <- data.frame("Mean"=colMeans(meanQ_HiseqX1),
                           "sd"=apply(meanQ_HiseqX1,2, sd),
                           "lane_read"=as.factor(colnames(meanQ_HiseqX1)),
                           "Machine"=as.factor("HiSeqX 1")) %>%
  melt(measure.vars=c("Mean"))

meanQ_HiseqX2 <- ExtractBotLvlTimeseries(c("mean_q"),"HiSeqX 2",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  select(grep("mean", colnames(.)))

meanQ_HiseqX2 <- data.frame("Mean"=colMeans(meanQ_HiseqX2),
                           "sd"=apply(meanQ_HiseqX2,2, sd),
                           "lane_read"=as.factor(colnames(meanQ_HiseqX2)),
                           "Machine"=as.factor("HiSeqX 2")) %>%
  melt(measure.vars=c("Mean"))

meanQ_HiseqX3 <- ExtractBotLvlTimeseries(c("mean_q"),"HiSeqX 3",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  select(grep("mean", colnames(.)))

meanQ_HiseqX3 <- data.frame("Mean"=colMeans(meanQ_HiseqX3),
                           "sd"=apply(meanQ_HiseqX3,2, sd),
                           "lane_read"=as.factor(colnames(meanQ_HiseqX3)),
                           "Machine"=as.factor("HiSeqX 3")) %>%
  melt(measure.vars=c("Mean"))


p1 <- ggplot(rbind(meanQ_HiseqX1, meanQ_HiseqX2, meanQ_HiseqX3),
       aes(y=value, x=lane_read, color=Machine)) +
  geom_point(aes(fill="Mean"), size=2.5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymax=value+sd, ymin=value-sd, width=0.2, fill="Std."), alpha=0.5, position=position_dodge(width=0.5)) +
  coord_flip() +
  theme_bw() +
  scale_color_brewer(type = "qual", palette=6) +
  scale_fill_manual("",
                    values=c("black","red")) +
  guides(fill=guide_legend(override.aes=list(shape=c(16,NA),
                                             linetype=c(0,12)))) +
  ylab("") +
  theme(legend.position="bottom")

##########################################################################################################################

###########################################################
Cycl <- c(140,160)
meanQ_HiseqX4 <- ExtractBotLvlTimeseries(c("mean_q"),"HiSeqX 4",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  select(grep("mean", colnames(.)))

meanQ_HiseqX4 <- data.frame("Mean"=colMeans(meanQ_HiseqX4),
                           "sd"=apply(meanQ_HiseqX4,2, sd),
                           "lane_read"=as.factor(colnames(meanQ_HiseqX4)),
                           "Machine"=as.factor("HiSeqX 4")) %>%
  melt(measure.vars=c("Mean"))

meanQ_HiseqX5 <- ExtractBotLvlTimeseries(c("mean_q"),"HiSeqX 5",All.df.reduced) %>%
  filter(cycles>=Cycl[1],
         cycles<=Cycl[2]) %>%
  # now remove those corresponding to rapid runs on the same cycle level, these only use 2 lanes and have NA on the rest
  na.omit() %>%
  select(grep("mean", colnames(.)))

meanQ_HiseqX5 <- data.frame("Mean"=colMeans(meanQ_HiseqX5),
                           "sd"=apply(meanQ_HiseqX5,2, sd),
                           "lane_read"=as.factor(colnames(meanQ_HiseqX5)),
                           "Machine"=as.factor("HiSeqX 5")) %>%
  melt(measure.vars=c("Mean"))

p2 <- ggplot(rbind(meanQ_HiseqX1, meanQ_HiseqX2,meanQ_HiseqX3,meanQ_HiseqX4,meanQ_HiseqX5),
             aes(y=value, x=lane_read, color=Machine)) +
  geom_point(aes(fill="Mean"), size=2.5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymax=value+sd, ymin=value-sd, width=0.2, fill="Std."), alpha=0.5, position=position_dodge(width=0.5)) +
  coord_flip() +
  theme_bw() +
  scale_color_brewer(type = "qual", palette=6) +
  scale_fill_manual("",
                    values=c("black","red")) +
  guides(fill=guide_legend(override.aes=list(shape=c(16,NA),
                                             linetype=c(0,12)))) +
  ylab("") +
  theme(legend.position="bottom")

grid.arrange(p1,p2,ncol=2)
