library(quantsel)
library(ggplot2)
library(dplyr)
library(parallel)

source("setup_demography.R")
CORES=4
reps=2

phenfile="phenos.csv"
afrqfile="afrq.csv"
dbfile="results1.db"


multithread=T
source("different_selection.R")

#This is how I do treatments when possible, this makes a fully factorial experiment
treats <- expand.grid(seedmix=c(0,0.1,0.5,1),
                      repro=c(0,0.25,0.5),
                      denstol=c(0,0.10,0.20),
                      K=5000,
                      reps=1:reps)

#treats <- treats[sample(nrow(treats),32),]
treats <- treats[sample(nrow(treats),8),]

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
               gens=100,sampint=20,rp=i,plt=FALSE)
    })
    phenos <- do.call(rbind,lapply(allreslst.eq,function(x)
    {
        as.data.frame(x[["phen"]])
    }))
    afreqs <- do.call(rbind,lapply(allreslst.eq,function(x)
    {
        as.data.frame(x[["afrq"]])
    }))
    Nn <- do.call(rbind,lapply(allreslst.eq,function(x)
    {
        as.data.frame(x[["Nn"]])
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
    
    if (!is.null(dbfile))
    {
        if ((j==0) & file.exists(dbfile)) file.remove(dbfile)
        con <- RSQLite::dbConnect(RSQLite::SQLite(),dbfile)
        res <- RSQLite::dbWriteTable(con,"phenos",phenos,append=T)
        res <- RSQLite::dbWriteTable(con,"afreqs",afreqs,append=T)
        res <- RSQLite::dbWriteTable(con,"Nn",Nn,append=T)
        RSQLite::dbDisconnect(con)
    }
    
}

print("finished simulation")

