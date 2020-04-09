"landscape.theta.k" <-
function(rland)
  {
    retval <- matrix(NA,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        if (dim(rland.tmp$individuals)[1]>1)
          for (j in 1:length(rland$loci))
            {
#              print(paste(i,j))
              alleledist <- as.factor(landscape.locus(lnum=j,rland.tmp)[,c(-1:-landscape.democol())])
              if (length(unique(alleledist))>1)
                retval[i,j] <- theta.k(alleledist)
              else
                retval[i,j] <- NA
            }
      }
    retval
  }

