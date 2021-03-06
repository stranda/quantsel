library(quantsel)
library(ggplot2)
library(dplyr)
source("setup_demography.R")



onerep <- function(plt=T)
{
    rland <- spartina.landscape(nloc=16, nphen=4)  #simulate 16 loci and 4 phenotypes
    expmat <- matrix(c(  #16rows for 16 loci, 4 cols for 4 phenotypes
        1,0,0,0, 
        1,0,0,0,
        1,0,0,0,
        1,0,0,0,
        0,1,0,0,
        0,1,0,0,
        0,1,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,1,0,
        0,0,1,0,
        0,0,1,0,
        0,0,0,1,
        0,0,0,1,
        0,0,0,1,
        0,0,0,1
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
                                          0,   0,   0, 0,   #long scale 
                                          0,   0,   0, 0,   #long shape  #no selection
                                          0,   0,   0, 0,   #mixture   #phenotype 2   #no selection
                                          0,   0,   0, 0,   #not used  #no selection
                                          0,   0,   0, 0,   #survive  
                                          0,   0,   0, 0,   #reproduce 
                                          0,   0,   0, 0    #density tolerance  #no selection
                                          ),
                                        ncol=4,byrow=T)
                                 )
    
    rland <- landscape.new.plasticity(rland)
    
    pl <- rland$plasticity
    pl[,1] <- pl[,1] * 0.7 #make everything  less that max for phenotype 1
    pl[1:10,1] <- 1 #set strips 1 through 10 to max

    rland <- landscape.new.plasticity(rland,pm=pl)
    
    rland <- landscape.new.phenohab(rland) #set up a no selection phenotype x habitat map
    
    
    initpopsize <- 250
    inits <- matrix(initpopsize,ncol=rland$intparam$habitats,nrow=2)
    
    rland <- landscape.new.individuals(rland,c(inits))
    
    l=landscape.simulate(rland,1) 
    
    if (plt==T) landscape.plot.phenotypes(l,1,F)
    
    
    phens=c(1,2,3,4) #represented as 0 in c++
    gen=15 #only really need one to show pattern
    sumlst=list()[1:gen]
    
                                        #pdf(paste0("gaps_",gapprop,".pdf"), width=15,height=7.5)
    
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

        save(file=paste0("plastic_int_res.rda"),sumlst,sumdf)
    }
