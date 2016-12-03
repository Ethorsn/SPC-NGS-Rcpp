library(mvnfast)
library(dplyr)
load("../Data/meanH.Rdata")
load("../Data/ICdata.Rdata")
load("../Data/Case1Mean.Rdata")
load("../Data/Case2Mean.Rdata")
Rcpp::sourceCpp("../FunctionsAndRcpp/RcppHotellingT2.cpp")
Rcpp::sourceCpp("../FunctionsAndRcpp/RcppMCUSUM.cpp")

M <- nrow(TransformedData)
p <- ncol(TransformedData)
alpha <- 0.01

QuantileHot <- p*(M-1)*(M+1)/((M-p)*M) * qf(1-alpha, p, M-p)

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M
#########################
ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(x,0,0),2),rep(0,length(mu0)-6))
}
ConstructMeanChartingStatistic <- function(data,k) {
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
h_mean_3 <- H_Listmean[[1]]$Intervals[nrow(H_Listmean[[1]]$Intervals),] %>% mean()
h_mean_4 <- H_Listmean[[2]]$Intervals[11,] %>% mean()
h_mean_5 <- H_Listmean[[3]]$Intervals[11,] %>% mean()

number <- 1e5
sc1 <- seq(0.2,1.5, by=0.05)

hot_res <- matrix(0,ncol=1, nrow = length(sc1))
cusm_res <- matrix(0,ncol=3,nrow = length(sc1))

for (j in 1:length(sc1)){
  for (i in 1:number){
    d_tmp <- rbind(mvnfast::rmvn(20,
                                 mu=mu0,
                                 sigma = Sigma0),
                   mvnfast::rmvn(1,
                                 mu=ChangeMeanFun(mu0,sc1[j]),
                                 sigma = Sigma0))
    
    h_tmp <- apply(d_tmp, 1, function(x) HotellingT2(x,mu0=mu0,Sigma0=Sigma0))
    cusm_tmp_3 <- ConstructMeanChartingStatistic(d_tmp,k=0.3)
    cusm_tmp_4 <- ConstructMeanChartingStatistic(d_tmp,k=0.4)
    cusm_tmp_5 <- ConstructMeanChartingStatistic(d_tmp,k=0.5)
    
    hot_res[j,1] <- hot_res[j,1] + ifelse(h_tmp[21]>QuantileHot,1,0)
    
    cusm_res[j,1] <- cusm_res[j,1] + ifelse(cusm_tmp_3[21]>h_mean_3, 1, 0)
    cusm_res[j,2] <- cusm_res[j,2] + ifelse(cusm_tmp_4[21]>h_mean_4, 1, 0)
    cusm_res[j,3] <- cusm_res[j,3] + ifelse(cusm_tmp_5[21]>h_mean_5, 1, 0)
  }
}

p_hotelling_sc1 <- (hot_res)/number
p_cusum_sc1 <- (cusm_res)/number

save(sc1, p_hotelling_sc1, p_cusum_sc1, file = "../Data/transient_change_sim.Rdata")
########## Scenario 2 Persistent change.
ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(-x,-x,x),2),rep(0,length(mu0)-6))
}

HotellingChange<- function(i){
  data <- MASS::mvrnorm(n=500, mu=ChangeMeanFun(mu0 = mu0, x=i), Sigma=Sigma0)
  Hotellings <- apply(data,1, function(j) HotellingT2(x_new=j, mu0=mu0, Sigma0=Sigma0))
  FirstAbove <- which(Hotellings>QuantileHot) %>% first()
  
  if (is.na(FirstAbove)==TRUE){
    FirstAbove <- 500
  }
  return(FirstAbove)
}
# Perform simulations for a limited set of the simulations values in the change vecotr.
exChange<- changeVector.case1[c(10,20,30,40,50)] 
Hot_list1 <- c()
for (i in exChange){
   mcuh <- replicate(10, HotellingChange(i=i))
   Hot_list1 <- c(Hot_list1, mean(mcuh))
}
save(Hot_list1, file="../Data/Scenario1Hotelling.Rdata")

ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(0,0,x),2),rep(0,length(mu0)-6))
}

HotellingChange<- function(i){
  data <- MASS::mvrnorm(n=500, mu=ChangeMeanFun(mu0 = mu0, x=i), Sigma=Sigma0)
  Hotellings <- apply(data,1, function(j) HotellingT2(x_new=j, mu0=mu0, Sigma0=Sigma0))
  FirstAbove <- which(Hotellings>QuantileHot) %>% first()
  
  if (is.na(FirstAbove)==TRUE)
  {
    FirstAbove <- 500
  }
  return(FirstAbove)
}

exChange<- changeVector.case2[c(1,3,5,8,10)] 
Hot_list2 <- c()
for (i in exChange){
  mcuh <- replicate(1e5, HotellingChange(i=i))
  Hot_list2 <- c(Hot_list2, mean(mcuh))
}
save(Hot_list2, file="../Data/Scenario2Hotelling.Rdata")