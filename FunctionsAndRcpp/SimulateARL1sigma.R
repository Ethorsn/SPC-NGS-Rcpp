### Simulate out of control behaviour
library(Rcpp)
library(dplyr)
library(RcppArmadillo)
sourceCpp("../FunctionsAndRcpp/RcppSimulateARL1sigma.cpp") 
load("../Data/ICdata.Rdata")
load("../Data/WishartH.Rdata")

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

M <- nrow(TransformedData)
mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M

Control.Limits.cov <- lapply(H_ListSigma, function(x) x$Intervals %>%
                                na.omit() %>%
                                tail(1) %>%
                                mean()) %>%
  unlist()

MYSEQ <- c(0.3,0.4,0.5)
N<- 1e5
##################################################### 
# Case 1 simulations
ChangeSigmaFun<- function(Sigma0,c) { 
  Sigma1 <- Sigma0
  diag(Sigma1)[1:6] <- diag(Sigma1)[1:6]+c 
  return(Sigma1)
}

changeVectorSigma.case1 <- seq(0.01,0.5, length=5) 

ARL1Sigma.case1 <- list()
##### Simulate ARL1

for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVectorSigma.case1){
    print("Start!")
    tmp <- SimulateARL1Sigma(N,Control.Limits.cov[i], k=MYSEQ[i], mu0=mu0, mu1 = mu0, n0=0, Sigma0 = Sigma0, Sigma1 =ChangeSigmaFun(Sigma0,j), 7)  
    Return.df <- c(Return.df,mean(tmp))
  }
  
  ARL1Sigma.case1 <- c(ARL1Sigma.case1,as.data.frame(Return.df))
  names(ARL1Sigma.case1)[i] <- MYSEQ[i]
  print("New k!")
}
ARL1Sigma.case1
EDSigma.case1 <- list()
##### Simulate ED
q <- 20
for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVectorSigma.case1){
    print("Start!")
    
    tmp <- SimulateARL1Sigma(N,Control.Limits.cov[i], k=MYSEQ[i], mu0=mu0, mu1 = mu0, n0=q, Sigma0 = Sigma0, Sigma1 =ChangeSigmaFun(Sigma0,j), 7)
    Return.df <- c(Return.df,mean(tmp))
  }
  
  EDSigma.case1 <- c(EDSigma.case1,as.data.frame(Return.df))
  names(EDSigma.case1)[i] <- MYSEQ[i]
  print("New k!")
}

save(ARL1Sigma.case1,EDSigma.case1,changeVectorSigma.case1, file="../Data/Case1Sigma.Rdata")

##################################################### 
# Case 2 simulations
ChangeSigmaFun<- function(Sigma0,c) { 
  Sigma1 <- Sigma0
  Sigma1[3,3] <- Sigma1[3,3]+c 
  Sigma1[6,6] <- Sigma1[6,6]+c 
  return(Sigma1)
}
changeVectorSigma.case2 <- seq(0.05,2.5, length=5) 

ARL1Sigma.case2 <- list()
##### Simulate ARL1

for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVectorSigma.case2){
    tmp <- SimulateARL1Sigma(N,Control.Limits.cov[i], k=MYSEQ[i], mu0=mu0, mu1 = mu0, n0=0, Sigma0 = Sigma0, Sigma1 = ChangeSigmaFun(Sigma0,j), 7)
    
    Return.df <- c(Return.df,mean(tmp))
  }
  Return.df <- as.data.frame(Return.df)
  ARL1Sigma.case2 <- c(ARL1Sigma.case2,Return.df)
  names(ARL1Sigma.case2)[i] <- MYSEQ[i]
  print("New k!")
}

EDSigma.case2 <- list()
##### Simulate ED
q <- 20
for (i in 1:(length(MYSEQ))){
  Return.df <- c()
  for (j in changeVectorSigma.case2){
    
    tmp <- SimulateARL1Sigma(N,Control.Limits.cov[i], k=MYSEQ[i], mu0=mu0, mu1 =mu0, n0=q, Sigma0 = Sigma0,Sigma1 = ChangeSigmaFun(Sigma0,j), 7)  
    Return.df <- c(Return.df,mean(tmp))
  }
  Return.df<- as.data.frame(Return.df)
  EDSigma.case2 <- c(EDSigma.case2,Return.df)
  names(EDSigma.case2)[i] <- MYSEQ[i]
  print("New k!")
}

save(ARL1Sigma.case2,EDSigma.case2,changeVectorSigma.case2, file="../Data/Case2Sigma.Rdata")
