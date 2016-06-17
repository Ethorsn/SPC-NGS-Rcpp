########################################################################################################################################
library(Rcpp)
library(RcppArmadillo)
# set working directory to source file locations or place in the 
sourceCpp('RcppSimulateARL0sigma.cpp')
load("../Data/ICdata.Rdata")

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

M <- nrow(TransformedData)
mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M

#################################################################################################################################
CalculateControlLimit <- function(k,mu0,Sigma0,interval_h, Nmax=30,M=5e5,ARL0=100,epsilon=1e-2, threads) {
  # Function which find the control limits h, for the CUSUM chart using the bisection algorithm and a absolute convergence criterion.
  #
  # Args:
  #    k: the allowance constant
  #    mu0: the in control mean
  #    Sigma0: the in control covariance matrix. It is assumed that this is positive definite
  #    interval_h : the interval to search for the control limits
  #    Nmax:  Maximum number of iterations to find h.
  #    M:  number of simulations to perform.
  #    ARL0: Target in control average run length 
  #    epsilon:  convergence criterion.
  #    threads: the number of threads to use in the simulations
  #
  # Returns:
  #    A list containing the simulated average run length and the interval considered in that iteration for each iteration
  
  intervals <- matrix(ncol=2,nrow=Nmax+1)
  intervals[1,] <- interval_h
  SimResult <- c()
  Cov <- Sigma0
  mu <- mu0
  for (i in 1:Nmax){
    # Create midpoints and replicate these for parrallel processing in CalculateARL0
    midpoint <- function(x,y) (x+y)/2
    h.start <- midpoint(interval_h[1],interval_h[2])
    # Simulate ARL0 for the h.start value
    Arls <- SimulateARL0Sigma(M,
                              h.start,
                              k,
                              mu,
                              Cov,
                              threads)
    
    ARLsim <- mean(Arls) 
    print("Simulated ARL0 equal to:")
    print(ARLsim)
    # absolute convergence criterion
    if(abs(ARLsim-ARL0)<epsilon){
      intervals[i+1,] <- rep(h.start,2)
      SimResult <- c(SimResult,ARLsim)
      break
    }else{
      if (ARLsim<ARL0) {
        # if we are below target
        interval_h[1] <- h.start
        print("Simulated ARL0 to small, increasing control limit.")
      }else{
      # if we are above target
        interval_h[2] <- h.start
        print("Simulated ARL0 to large, decreasing control limit.")
      }
    }
    intervals[i+1,] <- interval_h
    SimResult <- c(SimResult,ARLsim)
  }
  print("Simulations done")
  return(list(Intervals=intervals,
              ARLsims=SimResult))
}

H_ListSigma <- list()
MYSEQ <- c(0.3,0.4,0.5)
tmp <- list(CalculateControlLimit(k = MYSEQ[1],
                                  mu=mu0,
                                  Sigma0 = Sigma0,
                                  interval_h=c(0,2000),
                                  Nmax=20,
                                  M=1e5, 
                                  ARL0=100,
                                  epsilon=5e-2, 
                                  threads = 7))
print(tmp)
H_ListSigma <- c(H_ListSigma,tmp)
save(H_ListSigma,file="../Data/WishartH.Rdata")
###

tmp <- list(CalculateControlLimit(k = MYSEQ[2],
                                  mu=mu0,
                                  Sigma0 = Sigma0,
                                  interval_h=c(0,top2),
                                  Nmax=20,
                                  M=1e5, 
                                  ARL0=100,
                                  epsilon=5e-2, 
                                  threads = 7))

print(tmp)
H_ListSigma <- c(H_ListSigma,tmp)
save(H_ListSigma,file="../Data/WishartH.Rdata")
###
tmp <- list(CalculateControlLimit(k = MYSEQ[3],
                                  mu=mu0,
                                  Sigma0 = Sigma0,
                                  interval_h=c(0,top3),
                                  Nmax=20,
                                  M=1e5, 
                                  ARL0=100,
                                  epsilon=5e-2, 
                                  threads = 7))
print(tmp)
#load("../Data/WishartH.Rdata")
H_ListSigma <- c(H_ListSigma,tmp)
save(H_ListSigma,file="../Data/WishartH.Rdata")
