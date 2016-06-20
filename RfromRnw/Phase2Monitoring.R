## ----LoadLibsTransforms, chache=TRUE-------------------------------------
#Libraries
library(Rcpp)
library(RcppArmadillo)
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
# Load all data that is necessary
load("../Data/meanH.Rdata")
load("../Data/WishartH.Rdata")
load("../Data/ICdata.Rdata")
load("../Data/LambdaAndXtable.Rdata")
# IC parameters
K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

M <- nrow(TransformedData)

mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M

# R functions and scripts.
load("../Data/MySQL_data.Rdata")
All.df.reduced <- filter(All.df.reduced,year(Date)>2012)
source("../FunctionsAndRcpp/BotlvlSeriesExtraction.R")
# C++ Functions
# Rcpp::sourceCpp("../FunctionsAndRcpp/RcppCPDestimation.cpp")
# Not needed here
Rcpp::sourceCpp("../FunctionsAndRcpp/RcppHotellingT2.cpp",showOutput=FALSE, verbose = FALSE)
Rcpp::sourceCpp("../FunctionsAndRcpp/RcppMCUSUM.cpp",showOutput=FALSE, verbose = FALSE)
Rcpp::sourceCpp("../FunctionsAndRcpp/RcppWishart.cpp",showOutput=FALSE, verbose = FALSE)
# What types of cycles
Cycl <- c(102,126)
# what variables to extract
vars <- c("mean_q","pct_q30","error_rate")
###### Run on different type of run, cycles are different!!!!!!
Hiseq3.read <- ExtractBotLvlTimeseries(variable = vars,"HiSeq 3", All.df.reduced) %>% na.omit()
###### Same type etc
Hiseq4.read <- ExtractBotLvlTimeseries(variable = vars,"HiSeq 4", All.df.reduced) %>%
  filter(cycles < Cycl[2],cycles > Cycl[1]) %>%
  na.omit()
###### Same type etc
Hiseq5.read <- ExtractBotLvlTimeseries(variable = vars,"HiSeq 5", All.df.reduced) %>% 
  filter(cycles < Cycl[2],cycles > Cycl[1]) %>%
  na.omit()
#### Remove unwanted stuff
rm(t2,t4,df.sample.results,All.df,All.df.reduced,Cycl)

get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## ----TransformData-------------------------------------------------------
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
# K needs to be specified globally.
K <- 100
Hiseq3.Trans <- TransformAndDivide(Hiseq3.read[,-c(1:3)])
Hiseq4.Trans <- TransformAndDivide(Hiseq4.read[,-c(1:3)])
Hiseq5.Trans <- TransformAndDivide(Hiseq5.read[,-c(1:3)])

## ----ControlLimits-------------------------------------------------------
# extract M and p, number of observations and variables.
M <- nrow(TransformedData)
p <- ncol(TransformedData)
# set type 1 error
alpha <- 0.01
# Calculate 
QuantileHot <- p*(M-1)*(M+1)/((M-p)*M) * qf(1-alpha, p, M-p)
### Hotelling for the three different machines.
HiSeq3.Hots <- data.frame("HiSeq3"=apply(Hiseq3.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0)),
                          "x"=1:nrow(Hiseq3.Trans),
                          "Alarm"=as.factor(ifelse(apply(Hiseq3.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))>QuantileHot,1,0)))
HiSeq4.Hots <- data.frame("HiSeq4"=apply(Hiseq4.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0)),
                          "x"=1:nrow(Hiseq4.Trans),
                          "Alarm"=as.factor(ifelse(apply(Hiseq4.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))>QuantileHot,1,0)))
HiSeq5.Hots <- data.frame("HiSeq5"=apply(Hiseq5.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0)),
                          "x"=1:nrow(Hiseq5.Trans),
                          "Alarm"=as.factor(ifelse(apply(Hiseq5.Trans, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))>QuantileHot,1,0)))


## ----HiSeqPhase2Hotelling, fig.cap="Hotellings T-square statistic for the HiSeq 3 (left), HiSeq 4 (middle) and HiSeq 5 machines with in control parameters based on HiSeq 6 data. The straight horizontal line represents the control limit.", fig.height=4----
p1 <- ggplot() +
  geom_point(aes(x=x,y=HiSeq3, color=Alarm), data=HiSeq3.Hots) +
  geom_hline(yintercept=QuantileHot, alpha=0.6) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("Alarm", values = c("#4e4eff","#ff4e4e"), breaks=c("0","1")) +
  ggtitle("HiSeq 3") 

p2 <- ggplot() +
  geom_point(aes(x=x,y=HiSeq4, color=Alarm), data=HiSeq4.Hots) +
  geom_hline(yintercept=QuantileHot) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("Alarm", values = c("#4e4eff","#ff4e4e"), breaks=c("0","1")) +
  ggtitle("HiSeq 4")+
  theme(legend.position="none")

p3 <- ggplot() +
  geom_point(aes(x=x,y=HiSeq5, color=Alarm), data=HiSeq5.Hots) +
  geom_hline(yintercept=QuantileHot) +
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

lambda4 <- t(mu0-colMeans(Hiseq4.Trans))%*%Sigma0%*%(mu0-colMeans(Hiseq4.Trans))
lambda5 <- t(mu0-colMeans(Hiseq5.Trans))%*%Sigma0%*%(mu0-colMeans(Hiseq5.Trans))

M1 <- nrow(Hiseq4.Trans)
M2 <- nrow(Hiseq5.Trans)
detHiseq4<- det(var(Hiseq4.Trans)*(M1-1)/M1)
detHiseq5<- det(var(Hiseq5.Trans)*(M2-1)/M2)
detHiseq6<- det(Sigma0)

tmp <- data.frame("HiSeq 4" = detHiseq6/detHiseq4, "HiSeq 5" = detHiseq6/detHiseq5)
#xtable(tmp, label="tab:DetTable", caption="Table containing the Determinant of each estimated covaraince matrix")

## ----MCUSUM--------------------------------------------------------------
### MCUSUM
k <- 0.3
# Extract control limits for mean and covariance
covH <- H_ListSigma[[1]]$Intervals %>% 
  na.omit() %>% 
  tail(1) %>% 
  mean()

meanH <- H_Listmean[[1]]$Intervals %>%
  na.omit() %>%
  tail(1) %>%
  mean()
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
### Wishart transformation
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

#### MCUSUM mean
Mean_HiSeq4 <-  data.frame("HiSeq4"=ConstructMeanChartingStatistic(Hiseq4.Trans),
                           "x"=1:nrow(Hiseq4.Trans),
                           "Alarm"=as.factor(ifelse(ConstructMeanChartingStatistic(Hiseq4.Trans)>meanH,1,0)))

Mean_HiSeq5 <-  data.frame("HiSeq5"=ConstructMeanChartingStatistic(Hiseq5.Trans),
                           "x"=1:nrow(Hiseq5.Trans),
                           "Alarm"=as.factor(ifelse(ConstructMeanChartingStatistic(Hiseq5.Trans)>meanH,1,0)))
#### MCUSUM Wishart.
Wishart_Hiseq4 <- data.frame("HiSeq4"=ConstructWishartChartingStatistic(Hiseq4.Trans),
                             "x"=1:nrow(Hiseq4.Trans),
                             "Alarm"=as.factor(ifelse(ConstructWishartChartingStatistic(Hiseq4.Trans)>covH,1,0)))
Wishart_Hiseq5 <- data.frame("HiSeq5"=ConstructWishartChartingStatistic(Hiseq5.Trans),
                             "x"=1:nrow(Hiseq5.Trans),
                             "Alarm"=as.factor(ifelse(ConstructWishartChartingStatistic(Hiseq5.Trans)>covH,1,0)))

## ----HiSeq45MCUSUMfig, fig.cap="MCUSUM control charts monitoring the mean vector and covariance matrix. The in control parameters is based on transformed data from the HiSeq 6 machine. The horizontal line is the control limit calculated with k=0.3."----
p1 <- ggplot(Mean_HiSeq4) +
  geom_point(aes(x=x,y=HiSeq4, color=Alarm)) +
  geom_hline(yintercept=meanH) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual(" ", values = c("#4e4eff","#ff4e4e"), breaks=c("Alarm","No Alarm")) +
  ggtitle("HiSeq 4 - Mean")

p11 <- ggplot(Mean_HiSeq5) +
  geom_point(aes(x=x,y=HiSeq5, color=Alarm)) +
  geom_hline(yintercept=meanH) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual(" ", values = c("#4e4eff","#ff4e4e"), breaks=c("Alarm","No Alarm")) +
  ggtitle("HiSeq 5 - Mean")

p2 <- ggplot(Wishart_Hiseq4) +
  geom_point(aes(x=x,y=HiSeq4, color=Alarm)) +
  geom_hline(yintercept=covH) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("", values = c("#ff4e4e"), breaks=c("Alarm")) +
  ggtitle("HiSeq 4 - Covariance")  

p22 <- ggplot(Wishart_Hiseq5) +
  geom_point(aes(x=x,y=HiSeq5, color=Alarm)) +
  geom_hline(yintercept=covH) +
  ylab("value") +
  xlab("time") +
  theme_bw() +
  scale_color_manual("", values = c("#ff4e4e"), breaks=c("Alarm")) +
  ggtitle("HiSeq 5 - Covariance")  

grid.arrange(p1,p2, p11,p22, ncol=2)

