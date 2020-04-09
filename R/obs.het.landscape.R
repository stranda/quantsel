old.landscape.obs.het <- function (Rland) 
{
    tot <- 0
    rl <- matrix(NA, nrow = Rland$intparam$habitats, ncol = length(Rland$loci))
    for (j in unique(landscape.populations(Rland))) {
        for (loc in 1:length(Rland$loci)) {
 

            if (landscape.ploidy(Rland)[loc] == 1) {
                rl[j, loc] <- NA
            }
            else {
                freq.df <- data.frame(table(landscape.locus(loc, Rland)[landscape.populations(Rland) == 
                  j, c(-1:-9, -11)], landscape.locus(loc, Rland)[landscape.populations(Rland) == 
                  j, c(-1:-9, -10)]))

                rl[j, loc] <- (1 - sum(freq.df[as.character(freq.df[, 
                  1]) == as.character(freq.df[, 2]), 3])/sum(freq.df[, 
                  3]))
            }
        }
    }
    rl
}

landscape.obs.het <- function (Rland) 
{

    rl <- matrix(NA, nrow = Rland$intparam$habitats, ncol = length(Rland$loci))

    pl <- landscape.ploidy(Rland)
    strt <-  landscape.democol() + cumsum(pl) - pl + 1
    stp <- strt + (pl - 1)
    hetflag <- lapply(1:length(Rland$loci),function(x,ind,strt,stp,pops)
                      {
                        tmp <- cbind(aggregate(ind[,strt[x]]!=ind[,stp[x]],by=list(pop=pops),mean),rep(x,length(unique(pops))))
                        names(tmp) <- c("pop","het","loc")
                        tmp[,c(1,3,2)]
                      },
                      ind=Rland$individuals,
                      strt=strt,
                      stp=stp,
                      pops=landscape.populations(Rland)
                      )
    as.matrix(reshape(do.call(rbind,hetflag),direction="wide",timevar="loc",idvar="pop")[,-1])
}

