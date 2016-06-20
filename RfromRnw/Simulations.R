## ----LoadData------------------------------------------------------------
load("../Data/WishartH.Rdata")
load("../Data/meanH.Rdata")
load("../Data/Case1Mean.Rdata")
load("../Data/Case2Mean.Rdata")
load("../Data/Case1Sigma.Rdata")
load("../Data/Case2Sigma.Rdata")
load("../Data/ICdata.Rdata")
load("../Data/CPDsims.Rdata")
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

Rcpp::sourceCpp("../FunctionsAndRcpp/RcppHotellingT2.cpp")
###################### Found online function! webpage: http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization
get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
#######################

## ----ControlLimitAndTable, results='asis'--------------------------------
k <- c(0.3,0.4,0.5)
library(dplyr)
library(grid)
# Extract control limits for mean and covariance
Control.Limits.mean <- lapply(H_Listmean, function(x) x$Intervals %>%
  na.omit() %>%
  tail(1) %>%
  mean()) %>%
  unlist() %>%
  round(3)

Control.Limits.Cov <- lapply(H_ListSigma, function(x) x$Intervals %>%
  na.omit() %>%
  tail(1) %>%
  mean()) %>%
  unlist() %>%
  round(3)

tmp.xtab <- rbind(c("Mean",Control.Limits.mean), c("Covariance",Control.Limits.Cov)) %>%
  as.data.frame()
colnames(tmp.xtab) <- c("Type",k)

library(xtable)
#tmp.xtab.print <- xtable(tmp.xtab, digits = 2, caption="The control limits calculated using the function CalculateControlLimit for a set of allowance constants k.", label="Control")
#print.xtable(tmp.xtab.print,include.rownames=FALSE)

# The table have been slightly modified in appearance but NOT the content.

## ----HotellingsScenario1-------------------------------------------------

ChangeMeanFun<- function(mu0,x) { 
  mu0 + c(rep(c(-x,-x,x),2),rep(0,length(mu0)-6))
}
# 
HotellingChange<- function(i){
  data <- MASS::mvrnorm(n=500, mu=ChangeMeanFun(mu0 = mu0, x=i), Sigma=Sigma0)
  Hotellings <- apply(data,1, function(j) HotellingT2(x_new=j, mu0=mu0, Sigma0=Sigma0))
  FirstAbove <- which(Hotellings>QuantileHot) %>% first()
  
  
  if (is.na(FirstAbove)==TRUE)
  # did none detect the change? if yes, set to 500.
  {
    FirstAbove <- 500
  }
  return(FirstAbove)
}
# Perform simulations for a limited set of the simulations values in the change vecotr.
exChange<- changeVector.case1[c(10,20,30,40,50)] 
# Hot_list1 <- c()
# for (i in exChange){
#   mcuh <- replicate(1e5, HotellingChange(i=i))
#   Hot_list1 <- c(Hot_list1, mean(mcuh))
# }
# save(Hot_list1, file="../Data/Scenario1Hotelling.Rdata")
load("../Data/Scenario1Hotelling.Rdata")

## ----ARL1MeanCase1, fig.cap="Out-of-control ARL and ED for the MCUSUM chart of Scenario 1 - changes in the mean of all variables of lane one.", fig.height=4----
library(ggplot2)
library(gridExtra)
ggplot.tmp <- do.call("cbind",ARL1.case1) %>% as.data.frame()
ggplot.tmp$x <- changeVector.case1
colnames(ggplot.tmp) <- c(paste0("k",c(0.3,0.4,0.5)), "x")
p1 <- ggplot(ggplot.tmp)+
  geom_line(aes(y=k0.3,x=x, linetype="solid")) +
  geom_line(aes(y=k0.4,x=x, linetype="dashed")) +
  geom_line(aes(y=k0.5,x=x, linetype="dotdash")) +
  ylab(expression(ARL[1])) +
  xlab(expression(delta)) +
  scale_linetype("Allowance \n constant, k", labels=paste0("0.",3:5)) +
  theme_bw() +
  theme(legend.position="right") 
  

ggplot.tmp <- do.call("cbind",ED.case1) %>% as.data.frame()
ggplot.tmp$x <- changeVector.case1
colnames(ggplot.tmp) <- c(paste0("k",c(0.3,0.4,0.5)), "x")

p2 <- ggplot(ggplot.tmp) +
  geom_line(aes(y=k0.3,x=x), linetype=1) +
  geom_line(aes(y=k0.4,x=x), linetype=2) +
  geom_line(aes(y=k0.5,x=x), linetype=3) +
  theme_bw() +
  ylab("ED") +
  xlab(expression(delta)) +
  theme(legend.position="none") 
p3 <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

grid.arrange(p1,p2,p3, ncol=3, widths=c(2.3, 2.3, 0.8))

## ----StartConstructingXtbles---------------------------------------------
ARL1sigma <- do.call("rbind",ARL1Sigma.case1)
colnames(ARL1sigma) <-  changeVectorSigma.case1

EDsigma <- do.call("rbind", EDSigma.case1)
colnames(EDsigma) <-  changeVectorSigma.case1
EDsigma<- round(EDsigma,2)
#xtable(ARL1sigma)
#xtable(EDsigma)
# xtable code above generates the tables and was then pasted/manipulated into the printed format in the thesis. 

## ----HotellingT2case2----------------------------------------------------
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
#Hot_list2 <- c()
#for (i in exChange){
#  mcuh <- replicate(1e5, HotellingChange(i=i))
#  Hot_list2 <- c(Hot_list2, mean(mcuh))
#}
#save(Hot_list2, file="../Data/Scenario2Hotelling.Rdata")
load("../Data/Scenario2Hotelling.Rdata")

## ----ARL1meanCase2, fig.cap="Out-of-control ARL simulations of Scenario 2 - changes in the mean of the error rate of lane one.", fig.height=4----
ggplot.tmp <- do.call("cbind",ARL1.case2) %>% as.data.frame()
ggplot.tmp$x <- changeVector.case2
colnames(ggplot.tmp) <- c(paste0("k",c(0.3,0.4,0.5)), "x")
p1 <- ggplot(ggplot.tmp)+
  geom_line(aes(y=k0.3,x=x, linetype="solid")) +
  geom_line(aes(y=k0.4,x=x, linetype="dashed")) +
  geom_line(aes(y=k0.5,x=x, linetype="dotdash")) +
  ylab(expression(ARL[1])) +
  xlab(expression(delta)) +
  scale_linetype("Allowance \n constant", labels=paste0("0.",3:5)) +
  theme_bw() +
  theme(legend.position="right") 


ggplot.tmp <- do.call("cbind",ED.case2) %>% as.data.frame()
ggplot.tmp$x <- changeVector.case2
colnames(ggplot.tmp) <- c(paste0("k",c(0.3,0.4,0.5)), "x")

p2 <- ggplot(ggplot.tmp) +
  geom_line(aes(y=k0.3,x=x), linetype=1) +
  geom_line(aes(y=k0.4,x=x), linetype=2) +
  geom_line(aes(y=k0.5,x=x), linetype=3) +
  theme_bw() +
  ylab("ED") +
  xlab(expression(delta)) +
  theme(legend.position="none")
p3 <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

grid.arrange(p1,p2,p3, ncol=3, widths=c(2.3, 2.3, 0.8))

## ------------------------------------------------------------------------
ARL1sigma <- do.call("rbind",ARL1Sigma.case2)
colnames(ARL1sigma) <-  changeVectorSigma.case2

EDsigma <- do.call("rbind", EDSigma.case2)
colnames(EDsigma) <-  changeVectorSigma.case2
EDsigma<- round(EDsigma,2)

#xtable(ARL1sigma)
#xtable(EDsigma)

