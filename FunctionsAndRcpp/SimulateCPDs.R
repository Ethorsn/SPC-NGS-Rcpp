### Simulate out of control behaviour
library(Rcpp)
library(dplyr)
library(RcppArmadillo)
library(ggplot2)
# Sometimes you need to load another function with sourceCpp before loading this function,
# otherwise it does not recognize the compiler....
sourceCpp('RcppCPDestimation.cpp')

load("../Data/ICdata.Rdata")
load("../Data/Case1Mean.Rdata")
load("../Data/Case2Mean.Rdata")

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

M <- nrow(TransformedData)

mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M

##########################################
# Scenario 1
##########################################
ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(-x,-x,x),2),rep(0,length(mu0)-6))
}

# Take a few values of the sequence, not all.
CPDchangeVector.case1 <- seq(0.0001,0.05, length=50)[floor(seq(1,50,length=10))] 
# Repeat the process "number" times
number <- 1e5

Case1.CPD <- c()
# use ED k=0.3 to test CPD estimation.
for (j in 1:length(CPDchangeVector.case1)){
  print("Start!")
  tmp <- SimulateCPDs(NumbSims = number, N = floor(ED.case1$`0.3`[j]), tau = 20, mu0=mu0, mu1 = ChangeMeanFun(mu0, changeVector.case1[j]), sigma0 = Sigma0, sigma1 = Sigma0, No_threads=7)  
  Case1.CPD <- c(Case1.CPD,mean(tmp))
}
###########
# If covariance has changed. What is the result of D_bar?
###########
ChangeSigmaFun<- function(Sigma0,c) { 
  Sigma1 <- Sigma0
  diag(Sigma1)[1:6] <- diag(Sigma1)[1:6]+c 
  return(Sigma1)
}
ConstantSigma.case1 <- seq(0.01,0.5, length=5)[2]

Case1.CPD.ChangeSigma <- c()
for (j in 1:length(CPDchangeVector.case1)){
  print("Start!")
  tmp <- SimulateCPDs(NumbSims = number, N = floor(ED.case1$`0.3`[j]), tau = 20, mu0=mu0, mu1 = ChangeMeanFun(mu0, changeVector.case1[j]), sigma0 = Sigma0, sigma1 = ChangeSigmaFun(Sigma0,ConstantSigma.case2), No_threads=7)  
  Case1.CPD.ChangeSigma <- c(Case1.CPD.ChangeSigma,mean(tmp))
}
##########################################
# Scenario 2
##########################################
ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(0,0,x),2),rep(0,length(mu0)-6))
}

CPDchangeVector.case2 <- seq(0.05,2, length=10)
# Repeat the process "number" times

Case2.CPD <- c()
# use ED k=0.3 to test CPD estimation.
for (j in 1:length(CPDchangeVector.case2)){
  print("Start!")
  tmp <- SimulateCPDs(NumbSims = number, N = floor(ED.case2$`0.3`[j]), tau = 20, mu0=mu0, mu1 = ChangeMeanFun(mu0, CPDchangeVector.case2[j]), sigma0 = Sigma0, sigma1 = Sigma0, No_threads=7)  
  Case2.CPD <- c(Case2.CPD,mean(tmp))
}

ChangeSigmaFun<- function(Sigma0,c) { 
  Sigma1 <- Sigma0
  Sigma1[3,3] <- Sigma1[3,3]+c 
  Sigma1[6,6] <- Sigma1[6,6]+c 
  return(Sigma1)
}

ConstantSigma.case2 <- seq(0.05,2.5, length=5)[2]

Case2.CPD.ChangeSigma <- c()
for (j in 1:length(CPDchangeVector.case2)){
  print("Start!")
  tmp <- SimulateCPDs(NumbSims = number, N = floor(ED.case2$`0.3`[j]), tau = 20, mu0=mu0, mu1 = ChangeMeanFun(mu0, changeVector.case1[j]), sigma0 = Sigma0, sigma1 = ChangeSigmaFun(Sigma0,ConstantSigma.case2), No_threads=7)  
  Case2.CPD.ChangeSigma <- c(Case2.CPD.ChangeSigma,mean(tmp))
}



load("../Data/CPDsims.Rdata")
save(Case1.CPD, Case1.CPD.ChangeSigma, ConstantSigma.case1, CPDchangeVector.case1, Case2.CPD, Case2.CPD.ChangeSigma, ConstantSigma.case2, CPDchangeVector.case2, file="../Data/CPDsims.Rdata")


