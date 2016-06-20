library(knitr)
# 
purl("./EDA/ExploratoryDataAnalysis.Rnw", output="./RfromRnw/ExploratoryDataAnalysis.R")
purl("./Results/Phase2Monitoring.Rnw", output="./RfromRnw/Phase2Monitoring.R")
purl("./Results/Simulations.Rnw", output="./RfromRnw/Simulations.R")
