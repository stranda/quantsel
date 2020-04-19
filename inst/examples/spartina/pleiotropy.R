###
### the idea is to allow one phenotype affect fecundity and allow it to evolve the same
### way everywhere.  And one phenotype to affect density tolerance
### going to use phenotypes 1 and 2, respectively.  the others will drift
### each trait is determined by 2 genes that are independent and 2 that are in common (see expmat)
###


### interesting to note how important drift is in this example
### no plasticity

library(quantsel)
library(ggplot2)
library(dplyr)
source("setup_demography.R")

onerep <- function(plt=T)
{
    rland <- spartina.landscape(nloc=16, nphen=4, K=2000)  #simulate 16 loci and 4 phenotypes
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

    ph7=matrix(rep(c(1,2,0.4,1),rland$intparam$habitats),  #alpha=1, beta=2, range around 1 = 0.05 and use original sign (=0)
               nrow=rland$intparam$habitats,ncol=4,byrow=T)

    rland <- landscape.new.phenohab(rland,fitcomp=7,ph=ph7) #adding selection for phenotype 7 

     #note that the gradient goes the other way and there is less strength (range parameter)
    ph8=matrix(rep(c(1,2,0.3,0),rland$intparam$habitats),  #alpha=2, beta=1, range around 1 = 0.05 and use original sign (=0)
               nrow=rland$intparam$habitats,ncol=4,byrow=T) 

    rland <- landscape.new.phenohab(rland,fitcomp=8,ph=ph8) #adding selection for phenotype 7 (repro)

    initpopsize <- 250
    inits <- matrix(initpopsize,ncol=rland$intparam$habitats,nrow=2)
    #inits <- inits*0
    #inits[1,] <- 10
    rland <- landscape.new.individuals(rland,c(inits))
    
    l=landscape.simulate(rland,1) 
    
    if (plt==T) landscape.plot.phenotypes(l,1,F)
    
    
    phens=c(1,2,3,4) #represented as 0 in c++
    gen=50 #only really need one to show pattern
    sumlst=list()[1:gen]
    
    for (i in 1:gen)
    {
        print(dim(l$individuals))
        print(dim(l$individuals))
        if (dim(l$individuals)[1]>0) l.old=l
        l=landscape.simulate(l,1)
        if (((i %% 1)==0)&(plt==T))
        {
            par(mfrow=c(2,2))
            for (phen in phens)
                landscape.plot.phenotypes(l,phen,F)
            par(mfrow=c(1,1))
        }
        print(i)
        print(dim(l$individuals))
                                        #    print(landscape.allelefreq(l) )
        print(colMeans(landscape.phenotypes.c(l)))
        sumlst[[i]] <- list(phenosum=cbind(gen=i,data.frame(landscape.phenosummary(l))),
                            afreq=cbind(gen=i,landscape.allelefreq(l)))
            
        sumlst[[i]]$gen=i
    }
                                        #dev.off()
    sumlst
}

######################################## functions above execution below

if (!exists("multithread"))
    {

        lst <- onerep()  #run one replicate of the simulation, you could add parameters like num gens, dispersal, selection
                 #or just edit the function above

        sumdf <- do.call(rbind,lapply(lst,function(x) {x$phenosum}))

        save(file=paste0("two_trait.rda"),lst,sumdf)
    }
