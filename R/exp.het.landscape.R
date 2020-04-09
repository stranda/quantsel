
landscape.exp.het <- function (Rland) 
{
    pl <- landscape.ploidy(Rland)
    strt <-  9 + cumsum(pl) - pl + 1  
    stp <- strt + (pl - 1)
    rl <- matrix(0, nrow = Rland$intparam$habitats, ncol = length(Rland$loci))
    afrq <- landscape.allelefreq(Rland)
    homframe <- aggregate(afrq$Freq^2,by=list(pop=afrq$pop,loc=afrq$loc),sum)
    homframe$het <- 1-homframe$x
    homframe <- homframe[,-3]
    as.matrix(reshape(homframe,idvar="pop",timevar="loc",direction="wide")[,-1])
  }
