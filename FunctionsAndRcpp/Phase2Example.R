####################################################################
# This script aims to outline how the functions are to be run in which order and so on.
#
# What is presented here is a simplification and summarisation of what is presented in the thesis.
# The Figures will not look alike because of the change in control limits from the Phase1Example.R script.
# If one wants to replicate the Figures in the thesis it is suggested to source the Phase2Monitoring.R 
# file.
####################################################################

load("ExampleFromPhase1.Rdata")

# From loading  "ExampleFromPhase1.Rdata" we should have the following variables in our global enviorment:
######
# transformation parameters: Lambdas,
# In control mean and covariance: mu0,Sigma0,
# Control limits: h_cov,h_mean,h_Hot,
# Quality control variables: qq.variables
# Cycle setting: Cycl
# Instrument the control chart was created from: Instrument.variable
# The constant which was used in transformation: K
#####

source("../FunctionsAndRcpp/BotlvlSeriesExtraction.R")
# Load data by running this script 
# source("../FunctionsAndRcpp/mysqlScript.R")
# or load it directly if you do not have access to the database.
load("../Data/MySQL_data.Rdata")
####################################################################################
# Source functions which are going to be used in this Phase 2 example.
sourceCpp("../FunctionsAndRcpp/RcppCPDestimation.cpp",showOutput=FALSE, verbose = FALSE)
sourceCpp("../FunctionsAndRcpp/RcppHotellingT2.cpp",showOutput=FALSE, verbose = FALSE)
sourceCpp("../FunctionsAndRcpp/RcppMCUSUM.cpp",showOutput=FALSE, verbose = FALSE)
sourceCpp("../FunctionsAndRcpp/RcppWishart.cpp",showOutput=FALSE, verbose = FALSE)

####################################################################################
# What Instrument have we created the control charts for?
Instrument.variable

# What Cycle setting?
Cycl 

Hiseq3.read <- ExtractBotLvlTimeseries(variable = qq.variables,"HiSeq 3", All.df.reduced) %>% na.omit()
###### Same type etc
Hiseq4.read <- ExtractBotLvlTimeseries(variable = qq.variables,"HiSeq 4", All.df.reduced) %>%
  filter(cycles < Cycl[2],cycles > Cycl[1]) %>%
  na.omit()
###### Same type etc
Hiseq5.read <- ExtractBotLvlTimeseries(variable = qq.variables,"HiSeq 5", All.df.reduced) %>% 
  filter(cycles < Cycl[2],cycles > Cycl[1]) %>%
  na.omit()

#######################################################################################################
# Perform transformation procedure on new data using these two functions
TransformBoxCox <- function(x_new,lambda){
  # Function which transforms new data using a Box-Cox transformation, together with previously estimated lambdas.
  #
  # Args:
  #    x_new: A new observation or possibly observations. It is assumed that the observation is on the same order as the in control data
  #           If it is several, it needs to be on tidy format.
  #    lambda: A vector containing the transformation parameter for BoxCox transform.
  #
  # Returns: 
  #    The transformed x_new observation or matrix of observations using the Box-Cox transformation with parameter lambda
  
  # BoxCox transformation.
  if (nrow(x_new)>1) # is this one observation or a matrix of observations?
  {
    x_trans <- matrix(nrow=dim(x_new)[1], ncol=dim(x_new)[2]) 
    for (i in 1:ncol(x_new)){
      x_trans[,i]  <- unlist(forecast::BoxCox(x_new[,i],lambda = lambda[i]))
    }
  }
  else
  {
    x_trans <- c()
    for (i in 1:length(x_new)){
      x_trans[i]  <- forecast::BoxCox(x_new[i],lambda = lambda[i])
    }
  }
  return(x_trans)
}

TransformAndDivide<- function(x) {
  # Function which performs a box cox transformation for variables which does not pct in their names. 
  # For those variables which have pct in their name a transformation using the quantile function is performed.
  # 
  # Args: 
  #     x: A vector or possibly matrix/data frame containing observations on the same order as the in control sample.
  #        x needs to have column names and it needs to be equal to those used in the in control sample.
  # 
  # Returns: 
  #    a vector or possibly the matrix of transformed data. 
  
  tmp.names <- colnames(x)
  PCT <- grep("pct", tmp.names)
  
  tmp.ret <- x
  
  tmp.ret[,!1:ncol(x) %in% PCT] <- TransformBoxCox(x[,!1:ncol(x) %in% PCT], lambda = Lambdas)
  
  tmp.ret[,grep("mean",tmp.names)] <- tmp.ret[,grep("mean",tmp.names)]/K
  
  tmp.ret[,1:ncol(x) %in% PCT] <- apply(x[,1:ncol(x) %in% PCT],2,function(x) qnorm(x))
  
  return(tmp.ret)
}

Hiseq3.Trans <- TransformAndDivide(Hiseq3.read[,-c(1:3)])
Hiseq4.Trans <- TransformAndDivide(Hiseq4.read[,-c(1:3)])
Hiseq5.Trans <- TransformAndDivide(Hiseq5.read[,-c(1:3)])
### Hotelling for the three different machines.
HiSeq3.Hots <- data.frame("HiSeq3"=apply(Hiseq3.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0)),
                          "x"=1:nrow(Hiseq3.Trans),
                          "Alarm"=as.factor(ifelse(apply(Hiseq3.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))>h_Hot,1,0)))
HiSeq4.Hots <- data.frame("HiSeq4"=apply(Hiseq4.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0)),
                          "x"=1:nrow(Hiseq4.Trans),
                          "Alarm"=as.factor(ifelse(apply(Hiseq4.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))>h_Hot,1,0)))
HiSeq5.Hots <- data.frame("HiSeq5"=apply(Hiseq5.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0)),
                          "x"=1:nrow(Hiseq5.Trans),
                          "Alarm"=as.factor(ifelse(apply(Hiseq5.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))>h_Hot,1,0)))

##################################
# Figure in Thesis
library(ggplot2)
library(reshape2)
library(GGally)
library(dplyr)
library(gridExtra)
library(grid)
library(ggExtra)

get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 <- ggplot() +
  geom_point(aes(x=x,y=HiSeq3, color=Alarm), data=HiSeq3.Hots) +
  geom_hline(yintercept=h_Hot, alpha=0.6) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("Alarm", values = c("#4e4eff","#ff4e4e"), breaks=c("0","1")) +
  ggtitle("HiSeq 3") 

p2 <- ggplot() +
  geom_point(aes(x=x,y=HiSeq4, color=Alarm), data=HiSeq4.Hots) +
  geom_hline(yintercept=h_Hot) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("Alarm", values = c("#4e4eff","#ff4e4e"), breaks=c("0","1")) +
  ggtitle("HiSeq 4")+
  theme(legend.position="none")

p3 <- ggplot() +
  geom_point(aes(x=x,y=HiSeq5, color=Alarm), data=HiSeq5.Hots) +
  geom_hline(yintercept=h_Hot) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("Alarm", values = c("#4e4eff","#ff4e4e"), breaks=c("0","1")) +
  ggtitle("HiSeq 5") +
  theme(legend.position="none")

p4 <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")
grid.arrange(p1,p2,p3,p4,ncol=4, widths=c(2.3, 2.3, 2.3, 0.8))
rm(p1,p2,p3)

##########################################################################################
# Construct the charting statistics for the mean 

ConstructMeanChartingStatistic <- function(data) {
  # Function which calculates the charting statistic for the mean of Croisers MCUSUM control chart. 
  # Note that the following variables are set globally: mu0, Sigma0, k. It is assumed that the function MCUSUM is loaded.
  # 
  # Args:
  #    data: contains transformed data on the same format as the in control sample. That implies the same order of the columns.
  #
  # Returns:
  #    A vector of charting statistics for the mean.
  ret <- c()
  Sold <- rep(0, ncol(Sigma0))
  
  
  for(i in 1:nrow(data)){
    tmp <- MCUSUM(unlist(data[i,]),Sold = Sold, mu0=mu0, Sigma0 = Sigma0, k=k)
    
    Sold <- tmp$Snew
    Cstat <- tmp$Cstatistic
    ret <- c(ret,Cstat)
  }
  return(ret)
}
#### MCUSUM mean
Mean_HiSeq4 <-  data.frame("HiSeq4"=ConstructMeanChartingStatistic(Hiseq4.Trans),
                           "x"=1:nrow(Hiseq4.Trans),
                           "Alarm"=as.factor(ifelse(ConstructMeanChartingStatistic(Hiseq4.Trans)>h_mean,1,0)))

Mean_HiSeq5 <-  data.frame("HiSeq5"=ConstructMeanChartingStatistic(Hiseq5.Trans),
                           "x"=1:nrow(Hiseq5.Trans),
                           "Alarm"=as.factor(ifelse(ConstructMeanChartingStatistic(Hiseq5.Trans)>h_mean,1,0)))

##########################################################################################
# Construct the charting statistics for the covariance 

ConstructWishartChartingStatistic<- function(data) {
  # Function which construct the charting statistic for the covariance matrix using the properties of the singular Wishart distribution and Croisers MCUSUM.
  # Note that the following variables are set globally: mu0, Sigma0, k. It is also assumed that the following functions are loaded: MCUSUM, TransformObsWishart.
  #
  # Args: 
  #    data: A vector or matrix/data frame on the same format as the initial data which was used to estimate mu0 and sigma0.
  #
  # Returns: 
  #    A vector containing charting statitsics for the covariance matrix.
  library(dplyr)
  ret <- c()
  Sold <- rep(0, ncol(Sigma0)-1)
  
  for (j in 1:nrow(data)){
    # Perform the transformation for each column of the data (which is the same as the number of columns of Sigma0).
    tmp <- lapply(1:ncol(data),function(x) TransformObsWishart(observation=unlist(data[j,]),i=x,mu0=mu0,sigma0 = Sigma0))  
    # Take out the placement of the largest of the statistics
    WhichMaxStatistic <- lapply(tmp, function(x) MCUSUM(x,Sold=Sold,mu0=rep(0,ncol(Sigma0)-1),Sigma0 = diag(ncol(Sigma0)-1), k=k)$Cstatistic) %>% 
      unlist() %>% 
      which.max()
    # Perform the MCUSUM scheme again on the index which resulted in the largest value.
    tmp <- MCUSUM(tmp[[WhichMaxStatistic]],Sold=Sold ,mu0=rep(0,ncol(Sigma0)-1), Sigma0 = diag(ncol(Sigma0)-1), k=k)
    Sold <- tmp$Snew
    
    Hstatistic <- tmp$Cstatistic
    ret <- c(ret,Hstatistic)
  }
  return(ret)
}

#### MCUSUM Wishart.
Wishart_Hiseq4 <- data.frame("HiSeq4"=ConstructWishartChartingStatistic(Hiseq4.Trans),
                             "x"=1:nrow(Hiseq4.Trans),
                             "Alarm"=as.factor(ifelse(ConstructWishartChartingStatistic(Hiseq4.Trans)>h_cov,1,0)))
Wishart_Hiseq5 <- data.frame("HiSeq5"=ConstructWishartChartingStatistic(Hiseq5.Trans),
                             "x"=1:nrow(Hiseq5.Trans),
                             "Alarm"=as.factor(ifelse(ConstructWishartChartingStatistic(Hiseq5.Trans)>h_cov,1,0)))

#############################################
# Figures from the thesis. Note that it does not have the same control limits!!!!!
p1 <- ggplot(Mean_HiSeq4) +
  geom_point(aes(x=x,y=HiSeq4, color=Alarm)) +
  geom_hline(yintercept=h_mean) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual(" ", values = c("#4e4eff","#ff4e4e"), breaks=c("Alarm","No Alarm")) +
  ggtitle("HiSeq 4 - Mean")

p11 <- ggplot(Mean_HiSeq5) +
  geom_point(aes(x=x,y=HiSeq5, color=Alarm)) +
  geom_hline(yintercept=h_mean) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual(" ", values = c("#4e4eff","#ff4e4e"), breaks=c("Alarm","No Alarm")) +
  ggtitle("HiSeq 5 - Mean")

p2 <- ggplot(Wishart_Hiseq4) +
  geom_point(aes(x=x,y=HiSeq4, color=Alarm)) +
  geom_hline(yintercept=h_cov) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("", values = c("#ff4e4e"), breaks=c("Alarm")) +
  ggtitle("HiSeq 4 - Covariance")  

p22 <- ggplot(Wishart_Hiseq5) +
  geom_point(aes(x=x,y=HiSeq5, color=Alarm)) +
  geom_hline(yintercept=h_cov) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("", values = c("#ff4e4e"), breaks=c("Alarm")) +
  ggtitle("HiSeq 5 - Covariance")  

grid.arrange(p1,p2, p11,p22, ncol=2)

