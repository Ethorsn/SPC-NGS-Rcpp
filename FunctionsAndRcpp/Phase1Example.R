####################################################################
# This script aims to outline how the functions are to be run in which order and so on.
# It demands that OpenMP is enabled. For more information on how to install it, look in the README.
#
# What is presented here is a summarisation of what is performed in numerous scripts/.Rnw files.
# As done in the thesis we will focus HiSeq 6
####################################################################

# Run MySQL which fetches data and parses Instrument name and Date etc.
source("../FunctionsAndRcpp/mysqlScript.R")
# Remove 2012
df.sample.results <- filter(df.sample.results, year(Date)>2012)
All.df.reduced <- filter(All.df.reduced,year(Date)>2012)
######################################################################################
# This small section specifies the run settings we want to create the control charts for 

# What Instrument do we want to create the control chart for?
Instrument.variable <- "HiSeq 6"

# Cycle setting?
Cycl <- c(102,126)

# What variables are we going to construct control charts for?
qq.variables <- c("mean_q","pct_q30","error_rate")

# These settings do not include rapid or other "special runs"

#######################################################################################
# Extract in control data
Instrument.tag <- 
  filter(df.sample.results,
         Instrument == Instrument.variable, 
         cycles < Cycl[2],
         cycles > Cycl[1]) 

Instrument.read <- filter(All.df.reduced,
                      Instrument== Instrument.variable, 
                      cycles < Cycl[2],
                      cycles > Cycl[1])

CreateICsample <- function(x, y,quality_crit){
  # Function which creates a in control sample for tag and read level. 
  #
  # Args:
  #     x: A data frame which contains "tag level" measurements. 
  #        It is assumed that the variables specified in quality_crit are contained in this data frame.
  #     y: A data frame which contains "read level" measurements.
  #        It is assumed that the variables specified in quality_crit are contained in this data frame.
  #     quality_crit: A list containing the variables and limits which we would like to use as quality limits.
  #                   A example would be quality_crit=list(error_rate=1, pct_q30=0.7).
  #
  # Returns:
  #     Returns a list where each element is a data frame representing a in control data frame. 
  #     The list will have two element corresponding to tag and read level.
  library(dplyr)
  IC.read <- y %>%
    filter(error_rate < quality_crit$err,
           pf_clusters > quality_crit$pf_clusters) 
  
  # How many flowcells are removed?
  # FlowcellRemoved <- anti_join(IC.read,y,by="flowcell_id")$flowcell_id %>% unique() 
  
  # Remove flowcells which are not ok to the quality characteristics
  Data <- semi_join(x,IC.read,by="flowcell_id")
  
  return(list(IC.tag = Data,
              IC.read = IC.read))
}

IC.Instrument.sample <- CreateICsample(x=Instrument.tag, 
                                       y=Instrument.read,
                                       quality_crit = list(err=1, pf_clusters=1e8))
# source functions.
source("../FunctionsAndRcpp/BotlvlSeriesExtraction.R")
# Extract the data from this in control data set.
library(tidyr)
ICsample.Instrument <- ExtractBotLvlTimeseries(qq.variables,
                                               machine=Instrument.variable, 
                                               Data=group_by(IC.Instrument.sample$IC.read)) %>%
  ############################################################
  # Here we specify what type of run we would like to construct control charts for.
  # when using na.omit() we exclude runs which are not using all 8 lanes and both reads. 
  # In our example here we will only use v4 type runs for HiSeq machines.
  ############################################################
  na.omit()

#######################################################################################
# Transform in control data
TransformDataFun <- function(x){
  # Function which transforms data using a Box-Cox transformation. 
  # 
  # Args:
  #     x: a vector containing numeric data.
  #
  # Returns: 
  #     A list containing two elements. The first element of the list is the estimated transformation parameter. 
  #     The second is the transformed data.
  
  # The estimation of the transformation parameter.
  lambda <- forecast::BoxCox.lambda(x)
  
  return(list(Lambda=lambda,
              data=forecast::BoxCox(x,lambda)))
}


# First 3 columns are flowcell_id, cycles and date.
colnam <- colnames(ICsample.Instrument[,-c(1:3)])

# Any variable containing pct will be transformed using normal quantile function
PCT <- grep("pct",colnam)

# Prepare data frame.
TransformedData <- ICsample.Instrument[,-c(1:3)]
# Transform using Box-Cox
TransformedData.tmp <- apply(TransformedData[,!1:ncol(TransformedData) %in% PCT], 2, TransformDataFun)

Lambdas <- lapply(TransformedData.tmp,function(x) x$Lambda) %>% unlist()
# Extract Box-Cox transformaed data.
TransformedData.tmp <- lapply(TransformedData.tmp, function(x) x$data) %>% do.call("cbind",.)
# replace non-transformed data with transformed.
TransformedData[,!1:ncol(TransformedData) %in% PCT] <- TransformedData.tmp
# Transform pct variables with quantile function
TransformedData[,PCT] <- apply(TransformedData[,1:ncol(TransformedData) %in% PCT],2,function(x) qnorm(x))
####
#save(Lambdas, TransformedData, file="../Data/ICdata.Rdata")
####
rm(PCT,colnam,TransformedData.tmp)
#######################################################################################
# Test the transformed in control data. Does it show much evidence of being normally distributed?
#######################################################################################

K <- 100
OrgOrder <- colnames(TransformedData)
# divide by a constant!
TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)] <- TransformedData[,1:ncol(TransformedData) %in% grep("mean",OrgOrder)]/K

t2 <- MVN::hzTest(TransformedData)
t4 <- mvShapiroTest::mvShapiro.Test(as.matrix(TransformedData))
# Higher p-value is better!
#save(t2,t4,file = "../Data/LambdaAndXtable.Rdata")

#######################################################################################
# Estimate control limits
#######################################################################################

# Specify the target in control average run length, allowance constant, 
# number of simulations

#  This is just an example.

k <- 0.5 # ALLOWANCE CONSTANT
ARL0 <- 200 # TARGET
N <- 100 # NUMBER OF SIMS. SHOULD BE SET 10^6.
Nmax <- 4 #  NUMBER OF MAXIMUM STEPS
eps <- 5e-2 # CONVERGENCE CRITERION
INTERVAL_mean <- c(0,5000) # INTERVAL TO SEARCH FOR CONTROL LIMIT FOR THE MEAN!
INTERVAL_cov <- c(0,5)

alpha <- 0.01 # type 1 error

p <- ncol(TransformedData)
M <- nrow(TransformedData)

# in control mean and covariance
mu0 <- colMeans(TransformedData)
Sigma0 <- var(TransformedData)*(M-1)/M

################################################
# Find control limit for mean
library(Rcpp)
sourceCpp('../FunctionsAndRcpp/RcppSimulateARL0.cpp')

printf <- function(...) invisible(print(sprintf(...))) 

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
ControlLimit_mean <- CalculateControlLimit(k = k,
                                           mu=mu0,
                                           Sigma0 = Sigma0,
                                           interval_h=INTERVAL_mean,
                                           Nmax=Nmax,
                                           M=N, 
                                           ARL0=ARL0,
                                           epsilon=eps, 
                                           threads = 7)

h_mean <- ControlLimit_mean$Intervals %>%
  tail(1) %>%
  mean()

################################################
# Find control limit for covariance 
sourceCpp('RcppSimulateARL0sigma.cpp')

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

ControlLimit_cov <- CalculateControlLimit(k = k,
                                           mu=mu0,
                                           Sigma0 = Sigma0,
                                           interval_h=INTERVAL_cov,
                                           Nmax=Nmax,
                                           M=N, 
                                           ARL0=ARL0,
                                           epsilon=eps, 
                                           threads = 7)

h_cov <- ControlLimit_cov$Intervals %>%
  tail(1) %>%
  mean()

# Calculate Hotellings T^2 statistic 
h_Hot <- p*(M-1)*(M+1)/((M-p)*M) * qf(1-alpha, p, M-p)
#########
# The following needs to be saved.
save(Lambdas,mu0,Sigma0,h_cov,h_mean,h_Hot,qq.variables,Cycl,Instrument.variable,K, file="ExampleFromPhase1.Rdata")

# End. 

####################################################################
# Summarise:
# We have constructed a in control sample, estimated the in control parameters
# and then constructed control limits for Croisers MCUSUM and Hotellings T2 control charts.
####################################################################

