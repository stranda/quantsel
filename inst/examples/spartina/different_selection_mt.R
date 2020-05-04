library(quantsel)
library(ggplot2)
library(dplyr)
library(parallel)

source("setup_demography.R")
CORES=40
reps=10

phenfile="phenos.csv"
afrqfile="afrq.csv"
dbfile="results1.db"


multithread=T
source("different_selection.R")

#This is how I do treatments when possible, this makes a fully factorial experiment
treats <- expand.grid(seedmix=c(0,0.1,0.5,1),
                      repro=c(0,0.25,0.5),
                      denstol=c(0,0.25,0.5),
                      K=5000,
                      reps=1:reps)

#treats <- treats[sample(nrow(treats),32),]

##this line is important, onerep expects numeric inputs
for (i in 1:ncol(treats)) treats[,i] <- as.numeric(as.character(treats[,i]))

set.seed(Sys.time()) #any integer initializes the seed using time right now
treats$seed <- round(runif(nrow(treats),min=0,max=1000000000)) #these refer to RNG seeds
treats$subgrp <- 1:nrow(treats) %/% CORES

for (j in sort(unique(treats$subgrp)))
{
    print(j)
    allreslst.eq <- mclapply(which(treats$subgrp==j),mc.cores=CORES,function(i) {
        onerep(seedmix=treats$seedmix[i],
               repro=treats$repro[i],
               denstol=treats$denstol[i],
               K=treats$K[i],
               gens=20,rp=i,plt=FALSE)
    })
    phenos <- do.call(rbind,lapply(allreslst.eq,function(x)
    {
        as.data.frame(x[["phen"]])
    }))
    afreqs <- do.call(rbind,lapply(allreslst.eq,function(x)
    {
        as.data.frame(x[["afrq"]])
    }))

    for (m in 1:ncol(afreqs)) {afreqs[,m] <- as.numeric(as.character(afreqs[,m]))}
    
    if (j==0)
    {
        write.table(file=phenfile,sep=",",row.names=F,col.names=T,phenos)
        write.table(file=afrqfile,sep=",",row.names=F,col.names=T,afreqs)
    } else {
        write.table(file=phenfile,sep=",",row.names=F,col.names=F,append=T,phenos)
        write.table(file=afrqfile,sep=",",row.names=F,col.names=F,append=T,afreqs)
    }
}

print("finished simulation")

if (TRUE)
{
    tmp1 <- system(paste("rm -f",dbfile),intern=F)
    
    cmd <- paste0("csvsql --db sqlite:///",dbfile," --insert --tables phen ",phenfile)
    tmp1 <- system(paste0(cmd), intern=F)

    cmd <- paste0("csvsql --db sqlite:///",dbfile," --insert --tables afreq ",afrqfile)
    tmp1 <- system(paste0(cmd), intern=F)

}
