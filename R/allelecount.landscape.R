
landscape.genotypes <- function(Rland) #this function pastes the genotypes (haplotypes also) of all alleles
  {

    #make locus ids
    #

    
    tmp <- do.call(cbind,lapply(1:length(Rland$loci),function(x,ind,pl)
                                {
                                  strt <-  9 + cumsum(pl) - pl + 1
                                  stp <- strt + (pl - 1)
                                  if (pl[x]==1)
                                    {
                                      ind[,strt[x]]
                                    } else {
                                      paste(ind[,strt[x]],ind[,stp[x]],sep='/')
                                    }
                                },ind=Rland$individuals,pl=landscape.ploidy(Rland)
                                )
                   )
  }

"landscape.allelecount" <-
function(Rland,tbl.out=FALSE)
  {
    pl <- landscape.ploidy(Rland)
    strt <-  9 + cumsum(pl) - pl + 1
    stp <- strt + (pl - 1)

    counts <-    do.call(rbind,lapply(1:length(Rland$loci),function(i,Rland)
                                      {
                                        if (pl[i]==1)
                                          {
                                            tmp <- data.frame(table(pop=landscape.populations(Rland),
                                                                alleles=as.numeric(Rland$individuals[,strt[i]])))
                                          } else {
                                            tmp <- data.frame(table(pop=rep(landscape.populations(Rland),2),
                                                                    alleles=as.numeric(Rland$individuals[,strt[i]:stp[i]])))        
                                          }
                                        tmp <- cbind(tmp,loc=rep(i,dim(tmp)[1]))
                                      },Rland=Rland)
                         )

    rownames(counts) <- 1:dim(counts)[1]
    counts <- counts[,c(1,4,2,3)]
    if (tbl.out==T)
      {
        xtabs(Freq~pop+alleles+loc,counts)
      } else {
      counts
    }
  }

