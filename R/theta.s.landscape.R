"landscape.theta.s" <-
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
              retval[i,j] <- NA
              if ((rland$loci[[j]]$type==253))
                {
                  statevec <- landscape.locus.states(lnum=j,rland.tmp)$state
                  seqlen <- nchar(statevec[1])
                  segsites <- sum(apply(matrix(unlist(strsplit(statevec,split='')),ncol=seqlen,byrow=T),2,
                                        function (x){length(unique(x))})>1
                                  )
                  print(segsites)
                  if (segsites>0)
                    retval[i,j] <- theta.s(segsites,seqlen)[1]
                }
            }
      }
    retval
  }

