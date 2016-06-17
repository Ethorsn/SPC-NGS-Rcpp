#####################################################
# This script aims to find the control limit for Croisers MCUSUM for the mean!
#####################################################
library(Rcpp)
library(RcppArmadillo)
sourceCpp('../FunctionsAndRcpp/RcppSimulateARL0.cpp')
load("../Data/ICdata.Rdata")
printf <- function(...) invisible(print(sprintf(...))) 

OrgOrder <- colnames(TransformedData)
K <- 100 
TransformedData[,grep("mean",OrgOrder)] <- TransformedData[,grep("mean",OrgOrder)]/K

M <- nrow(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M
mu0 <- colMeans(TransformedData)

#################################################################################################################################
CalculateControlLimit <- function(k,mu0,Sigma0,interval_h, Nmax=1e5,M=1e6,ARL0=20,epsilon=1e-6, threads) {
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
  
  for (i in 1:Nmax){
    # Create midpoints and replicate these for parrallel processing in CalculateARL0
    midpoint <- function(x,y) (x+y)/2
    h.start <- midpoint(interval_h[1],interval_h[2])
    # Simulate ARL0 for the h.start value
    Arls <- SimulateARL0(M,
                         h.start,
                         k,
                         mu0,
                         Sigma0,
                         threads)
    
    ARLsim <- mean(Arls) 
    printf("Simulater ARL0 equal to: %f", ARLsim)
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

H_Listmean <- list()
MYSEQ <-c(0.3,0.4,0.5)
tmp <- list(CalculateControlLimit(k = MYSEQ[1],
                                  mu=mu0,
                                  Sigma0 = Sigma0,
                                  interval_h=c(0,5000),
                                  Nmax=40,
                                  M=5e5, 
                                  ARL0=100,
                                  epsilon=5e-2, 
                                  threads = 7))
H_Listmean  <- c(H_Listmean ,tmp)

tmp <- list(CalculateControlLimit(k = MYSEQ[2],
                                  mu=mu0,
                                  Sigma0 = Sigma0,
                                  interval_h=c(0,2580),
                                  Nmax=40,
                                  M=5e5, 
                                  ARL0=100,
                                  epsilon=5e-2, 
                                  threads = 7))
print(tmp)
H_Listmean <- c(H_Listmean ,tmp)
#save(H_Listmean, file="../Data/meanH.Rdata")
H_Listmean

#load("../Data/meanH.Rdata")
tmp <- list(CalculateControlLimit(k = MYSEQ[3],
                                  mu=mu0,
                                  Sigma0 = Sigma0,
                                  interval_h=c(0,2101),
                                  Nmax=40,
                                  M=5e5, 
                                  ARL0=100,
                                  epsilon=5e-2, 
                                  threads = 7))
print(tmp)
H_Listmean <- c(H_Listmean ,tmp)
#save(H_Listmean, file="../Data/meanH.Rdata")
