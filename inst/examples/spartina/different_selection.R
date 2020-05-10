###
### the idea is to allow one phenotype affect fecundity and allow it to evolve the same
### way everywhere.  And one phenotype to affect density tolerance
### going to use phenotypes 1 and 2, respectively.  the others will drift
### In this case there is no pleiotropy.
###
### Here selection is different in three zones
### 8 strips near the creek (pops 1-5)
### 8 strips intermediate
### 34 strips upland

### no plasticity

library(quantsel)
library(ggplot2)
library(dplyr)
source("setup_demography.R")

onerep <- function(seedmix=0.001,repro=0.4,denstol=0.1, K=1000, gens=100,rp=1, sampint=10, plt=T)
{
    rland <- spartina.landscape(nloc=16, nphen=4, K=K,seedmix)  #simulate 16 loci and 4 phenotypes
    expmat <- matrix(c(  #16rows for 16 loci, 4 cols for 4 phenotypes
        1  ,  0,0,0, 
        1  ,  0,0,0,
        0.5,  0.5,0,0,
        0.5,  0.5,0,0,
        0.5,  0.5,0,0,
        0.5,  0.5,0,0,
        0  ,    1,0,0,
        0  ,    1,0,0,
        0  ,    0,1,0,
        0  ,    0,1,0,
        0  ,    0,1,0,
        0  ,    0,1,0,
        0  ,    0,0,1,
        0  ,    0,0,1,
        0  ,    0,0,1,
        0  ,    0,0,1
    ),byrow=T,ncol=4)
    hsq <- c(1,1,1,1)
    rland <- landscape.new.expression(rland,
                                      expmat=expmat*0.125, #0.125 -> 1 diploid locus per phen, 
                                      hsq=hsq) #up to 8 alelle additive doses, when summed across 4 loci.
                                        #this standardizes the phenotype to range from 0-1
    rland <- landscape.new.gpmap(rland,
                                 ## 4 cols 5 rows.  Cols correspond to phenotype effects on fit components
                                 ##for each phenotype (0 is no effect, 4 phenotypes in this example)
                                 ##phenotypes are in C indexing so, add 1 to compare to pehnotypes above
                                 matrix(c(0,   0,   0, 0,   #short scale #no selection
                                          0,   0,   0, 0,   #long scale #no selection
                                          0,   0,   0, 0,   #long shape  #no selection
                                          0,   0,   0, 0,   #mixture   #phenotype 2   #no selection
                                          0,   0,   0, 0,   #not used  #no selection
                                          0,   0,   0, 0,   #survive  #no selection
                                          1,   0,   0, 0,   #reproduce 
                                          0,   1,   0, 0    #density tolerance  
                                          ),
                                        ncol=4,byrow=T)
                                 )
    
    rland <- landscape.new.plasticity(rland)
    
    pl <- rland$plasticity
    
    rland <- landscape.new.plasticity(rland,pm=pl)
    
    rland <- landscape.new.phenohab(rland) #set up a no selection phenotype x habitat map

    ph7=matrix(c(rep(c(1,2,repro,0),8),  #strip closest to creek
                 rep(c(1,1,0,0),8),    #intermediate
                 rep(c(2,1,repro,0),rland$intparam$habitats-16)), #rest of strips
               nrow=rland$intparam$habitats,ncol=4,byrow=T)

    rland <- landscape.new.phenohab(rland,fitcomp=7,ph=ph7) #adding selection for phenotype 7 

     #note that the gradient goes the other way in each location and there is less strength (range parameter)
    ph8=matrix(c(rep(c(2,1,denstol,0),8),  #strip closest to creek
                 rep(c(1,1,0,0),8),    #intermediate
                 rep(c(1,2,denstol,0),rland$intparam$habitats-16)), #rest of strips
               nrow=rland$intparam$habitats,ncol=4,byrow=T)

    rland <- landscape.new.phenohab(rland,fitcomp=8,ph=ph8) #adding selection for phenotype 8 (dens)

    initpopsize <- 400
    inits <- matrix(initpopsize,ncol=rland$intparam$habitats,nrow=2)
    #inits <- inits*0
    #inits[1,] <- 10
    rland <- landscape.new.individuals(rland,c(inits))
    
    l=landscape.simulate(rland,1) 
    
    if (plt==T) landscape.plot.phenotypes(l,1,F)
    
    
    phens=c(1,2,3,4) #represented as 0 in c++


    sgens <- (gens%/%sampint) 

    sumlst=list()[1:sgens]    
    for (i in 1:sgens)
    {
        print(dim(l$individuals))
        if (dim(l$individuals)[1]>0) l.old=l
        l=landscape.simulate(l,ifelse(i==sgens,gens%%sampint,sampint))
        if (plt==T)
        {
            par(mfrow=c(2,2))
            for (phen in phens)
                landscape.plot.phenotypes(l,phen,F)
            par(mfrow=c(1,1))
        }
        
        print(paste("Sample:",i,"Generation:",l$intparam$currentgen))
                                        #    print(landscape.allelefreq(l) )
        print(colMeans(landscape.phenotypes.c(l)))
        gn <- l$intparam$currentgen
        sumlst[[i]] <- list(phenosum=cbind(gen=gn,data.frame(landscape.phenosummary(l))),
                            afreq=cbind(gen=gn,landscape.allelefreq(l)),
                            neigh=landscape.neighborhood(l),
                            treat=cbind(seedmix=seedmix,repro=repro,denstol=denstol))
            
        sumlst[[i]]$gen=gn
    }
                                        #dev.off()
    phen <- do.call(rbind,lapply(sumlst,function(x){as.data.frame(x[["phenosum"]])}))
    phen$rep <- rp
    phen <- cbind(data.frame(sumlst[[1]]$treat[rep(1,nrow(phen)),]),phen)
    print("phen done")
    afrq <- do.call(rbind,lapply(sumlst,function(x){as.data.frame(x[["afreq"]])}))
    print(dim(afrq))
    afrq$rep <- rp
    afrq <- cbind(data.frame(sumlst[[1]]$treat[rep(1,nrow(afrq)),]),afrq)
    Nn <- do.call(rbind,lapply(sumlst,function(x){g=x[["gen"]];data.frame(gen=g,pop=1:length(x[["neigh"]]$Nn),Nn=x[["neigh"]]$Nn)}))
    list(phen=phen,afrq=afrq,Nn=Nn,landscape=l)
}

######################################## functions above execution below

if (!exists("multithread"))
    {

        lst <- onerep(seedmix=0.02,gens=500,sampint=50)  #run one replicate of the simulation, you could add parameters like num gens, dispersal, selection
                 #or just edit the function above

        sumdf <- do.call(rbind,lapply(lst,function(x) {x$phenosum}))

        save(file=paste0("two_trait.rda"),lst,sumdf)
    }
