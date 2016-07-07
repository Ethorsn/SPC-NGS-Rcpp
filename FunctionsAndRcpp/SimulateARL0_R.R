MCUSUM_R <- function(S_old,X_new,mu0,Sigma0,k){
  # Short function which calculates Croisers MCUSUM scheme.
  
  # Args:
  #   S_old: Vector from old iteration
  #   X_new: new observation
  #   mu0: in control mean vector
  #   Sigma0: in control covariance matrix.
  #   k: allowance constant
  S_old <- as.vector(S_old)
  X_new <- as.vector(X_new)
  mu0 <- as.vector(mu0)
  Sigma0 <- as.matrix(Sigma0)
  
  C <- t(S_old+X_new-mu0)%*%solve(Sigma0)%*%(S_old+X_new-mu0)
  
  if (C <= k){
    S_new <- rep(0,length(mu0))
  }else{
    S_new <- (S_old+X_new-mu0)*(1-k/C)
  }
  
  H <- t(S_new)%*%solve(Sigma0)%*%(S_new)
  
  return(list("H"=H,
              "S_old"=S_new))
}
# library(mvtnorm)
# test_obs <- rmvnorm(1,mean=mu0,sigma=Sigma0)
# MCUSUM_R(S_old=rep(0,length(mu0)), X_new=test_obs, mu0=mu0, Sigma0 = Sigma0, k=0.5)

SimulateArl0_R <- function(n, h, k, mu0, Sigma0, No_threads){
  
  # Load packages needed.
  library(foreach)
  library(doMC)
  library(mvtnorm)
  
  # Small function which simulates one RL 
  Sim_RL <- function(h=h,k=k,mu0=mu0,sigma0=Sigma0){
    # Large number to ensure that we dnt get stuck in a long loop
    LargeNumber <- 10^5
    # Inits.
    Count <- 1
    H <- 0
    S_old <- rep(0,length(mu0))
    # For loop to simulate run length
    for (i in 1:LargeNumber){
      X_new = rmvnorm(1,mean = mu0, sigma = Sigma0)
      ret <- MCUSUM_R(S_old = S_old, X_new = X_new, mu0 = mu0, Sigma0 = Sigma0, k=k)  
      
      if (ret$H>h) break
      
      S_old <- ret$S_old
      Count <- Count + 1
    }
    return(Count)
  }
  
  if (No_threads>1){
    # in parrallel
    registerDoMC(cores=No_threads)
    Return_vec <- times(n) %dopar% Sim_RL(h=h,k=k,mu0,Sigma0)
  }else{
    # in serial
    Return_vec <- replicate(n, Sim_RL(h=h,k=k,mu0,Sigma0))  
  }
  return(Return_vec)
}


