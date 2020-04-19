library(quantsel)
library(ggplot2)
library(dplyr)
library(parallel)

source("setup_demography.R")
CORES=16
reps=CORES*2

reps=2

multithread=T
source("different_selection.R")


#This is how I do treatments when possible, this makes a fully factorial experiment
treats <- expand.grid(seedmix=c(0,0.1,0.5,1),
                      repro=c(0,0.1,0.2,0.5),
                      denstol=c(0,0.1,0.2,0.5),
                      reps=1:reps)

set.seed(123123) #any integer initializes the seed
treats$seed <- round(runif(nrow(treats),min=0,max=1000000000)) #these refer to RNG seeds


allreslst.eq <- mclapply(1:nrow(treats),mc.cores=CORES,function(i) {
    print(i)
    set.seed(treats$seed[i])
    onerep(treats$seedmix[i],treats$repro[i],treats$denstol[i],treats$K[i],gens=50,rep=i,plt=FALSE)
})

print("finished simulation")
allres <- simsum(allreslst.eq,CORES=5) #gather up results into a list with 3 components

###tidyverse version of a quick visualization of phenotypes

plotdf <- allres$phen
plotdf$loc <- ifelse(plotdf$pop<9,"bank",ifelse(plotdf$pop>16,"rest","intermediate"))

pdf("diffselection.pdf")

for (r in unique(plotdf$repro))
    for (d in unique(plotdf$denstol))
        {
            df  <- plotdf[(plotdf$repro==r)&(plotdf$denstol==d),]
            print(
                ggplot(df,aes(x=gen,y=mean,color=phenotype)) + geom_smooth() +
                facet_wrap(~loc+seedmix) +
                ylim(c(0,1)) + ggtitle(paste0("selection on two traits when repro = ",r," and denstol = ",d))
                )
        }

dev.off()

