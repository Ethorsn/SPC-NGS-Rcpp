library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../FunctionsAndRcpp/RcppSimulateARL0.cpp")
load("../Data/meanH.Rdata")
load("../Data/ICdata.Rdata")
source("./SimulateARL0_R.R")

h <- H_Listmean[[1]]$Intervals %>% na.omit() %>% tail(1) %>% mean() %>% round(2)
k <- 0.3

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K
M <- nrow(TransformedData)
mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M

benchman <- microbenchmark(
  #1e2 simulations
  "1e2_Rcpp_OpenMP"=SimulateARL0(1e2, h,k,mu0,Sigma0,7),
  "1e2_R"=SimulateArl0_R(1e2, h, k ,mu0,Sigma0,No_threads = 1),
  "1e2_R_foreach"=SimulateArl0_R(1e2, h, k ,mu0,Sigma0, No_threads = 7),
  #1e3 simulations
  "1e3_Rcpp_OpenMP"=SimulateARL0(1e3, h,k,mu0,Sigma0,7),
  "1e3_R"=SimulateArl0_R(1e3, h,k,mu0,Sigma0,1),
  "1e3_R_foreach"=SimulateARL0(1e3, h,k,mu0,Sigma0,7),
  #1e4 simulations
  "1e4_Rcpp_OpenMP"=SimulateARL0(1e4, h,k,mu0,Sigma0,7),
  "1e4_R"=SimulateArl0_R(1e4, h,k,mu0,Sigma0,1),
  "1e4_R_foreach"=SimulateArl0_R(1e4, h,k,mu0,Sigma0,7),
  #1e5 simulations
  "1e5_Rcpp_OpenMP"=SimulateARL0(1e5, h,k,mu0,Sigma0,7),
  "1e5_R"=SimulateArl0_R(1e5, h,k,mu0,Sigma0,1),
  "1e5_R_foreach"=SimulateArl0_R(1e5, h,k,mu0,Sigma0,7),
  
  times=10)
save(benchman, file = "../Data/Benchmarks.Rdata")
autoplot(benchman)
