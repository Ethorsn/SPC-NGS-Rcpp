############################################################
# This script simulates out-of-control behaviour for the mean control chart.
# It is assumed that the control limits have been determined.
############################################################
library(Rcpp)
library(dplyr)
library(RcppArmadillo)
# Load function
sourceCpp("../FunctionsAndRcpp/RcppSimulateARL1.cpp") 
# load in control data and control limits for the mean.
load("../Data/ICdata.Rdata")
load("../Data/meanH.Rdata")

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)

Control.Limits.mean <- lapply(H_Listmean, function(x) x$Intervals %>%
                                na.omit() %>%
                                tail(1) %>%
                                mean()) %>%
  unlist()

MYSEQ <- c(0.3,0.4,0.5)
N<- 1e5
##################################################### 
# Scenario 2 simulations
ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(-x,-x,x),2),rep(0,length(mu0)-6))
}

changeVector.case1 <- seq(0.0001,0.05, length=50) 

ARL1.case1 <- list()
##### Simulate ARL1
for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVector.case1){
    h <- Control.Limits.mean[i]
    k <- MYSEQ[i]
    tmp <- SimulateARL1mean(n=N, h=h, k=k, mu0=mu0, mu1 = ChangeMeanFun(mu0, j), n0=0, Sigma0 = Sigma0, 7)  
    
    Return.df <- c(Return.df,mean(tmp))
  }
  
  ARL1.case1<- c(ARL1.case1,as.data.frame(Return.df))
  names(ARL1.case1)[i] <- MYSEQ[i]
  print("New k!")
}

ED.case1 <- list()
q <- 20
##### Simulate ED
for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVector.case1){
    tmp <- SimulateARL1mean(N,Control.Limits.mean[i], k=MYSEQ[i], mu0=mu0, mu1 = ChangeMeanFun(mu0, j), n0=q, Sigma0 = Sigma0, 7)  
    
    Return.df <- c(Return.df,mean(tmp))
  }
  
  ED.case1<- c(ED.case1,as.data.frame(Return.df))
  names(ED.case1)[i] <- MYSEQ[i]
  print("New k!")
}

save(ARL1.case1,ED.case1,changeVector.case1, file="../Data/Case1Mean.Rdata")

##################################################### 
# Scenario 3 simulations
ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(0,0,x),2),rep(0,length(mu0)-6))
}

changeVector.case2 <- seq(0.01,2, length=50) 

ARL1.case2 <- list()
##### Simulate ARL1

for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVector.case2){
    tmp <- SimulateARL1mean(N,Control.Limits.mean[i], k=MYSEQ[i], mu0=mu0, mu1 = ChangeMeanFun(mu0, j), n0=0, Sigma0 = Sigma0, 7)  
    
    Return.df <- c(Return.df,mean(tmp))
  }
  Return.df <- as.data.frame(Return.df)
  names(ARL1.case2)[i] <- MYSEQ[i]
  ARL1.case2 <- c(ARL1.case2,Return.df)
  print("New k!")
}

ED.case2 <- list()
q <- 20
##### Simulate ED
for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVector.case2){
    tmp <- SimulateARL1mean(N,Control.Limits.mean[i], k=MYSEQ[i], mu0=mu0, mu1 = ChangeMeanFun(mu0, j), n0=q, Sigma0 = Sigma0, 7)  
    
    Return.df <- c(Return.df,mean(tmp))
  }
  
  ED.case2 <- c(ED.case2,as.data.frame(Return.df))
  names(ED.case2)[i] <- MYSEQ[i]
  print("New k!")
}
save(ARL1.case2,ED.case2,changeVector.case2, file="../Data/Case2Mean.Rdata")
