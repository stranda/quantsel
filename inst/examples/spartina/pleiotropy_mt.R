library(quantsel)
library(ggplot2)
library(dplyr)
library(parallel)

source("setup_demography.R")
CORES=16
reps=CORES*2
multithread=T
source("pleiotropy.R")

set.seed(123124) #any integer initializes the seed
seeds <- round(runif(reps,min=0,max=1000000000)) #these refer to RNG seeds
allreslst.eq <- mclapply(1:reps,mc.cores=CORES,function(i) {print(i);set.seed(seeds[i]);onerep(plt=FALSE)[[1]]})
print("finished simulation")
allres <- simsum(allreslst.eq,CORES=CORES)
