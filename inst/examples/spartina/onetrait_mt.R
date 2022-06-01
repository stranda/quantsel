library(quantsel)
library(ggplot2)
library(dplyr)
library(parallel)

source("setup_demography.R")
CORES=16
reps=CORES*2
multithread=T
source("onetrait.R")


allreslst.eq <- mclapply(1:reps,mc.cores=CORES,function(i) {print(i);onerep(plt=FALSE)[[1]]})
print("finished simulation")
allres <- simsum(allreslst.eq,CORES=CORES)

