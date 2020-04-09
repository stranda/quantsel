"landscape.tajima.d" <-
function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            retval[i,j] <- NA
            if ((rland$loci[[j]]$type==253))
              {
                alleledist <- as.factor(landscape.locus(lnum=j,rland.tmp)[,c(-1,-2,-3)])
                if (length(unique(alleledist))>1)
                  theta.ewens <- c(theta.k(alleledist),0)
                else
                  theta.ewens <- c(NA,NA)
                statevec <- landscape.locus.states(lnum=j,rland.tmp)$state
                seqlen <- nchar(statevec[1])
                                        #            print(j)
                                        #                print(paste("len statevec",length(statevec)))
                                        #            print(seqlen)
                segsites <- sum(apply(matrix(unlist(strsplit(statevec,split='')),ncol=seqlen,byrow=T),2,function (x){length(unique(x))})>1)
                if (segsites>0)
                  retval[i,j] <- (theta.ewens[1]-theta.s(segsites,seqlen)[1])/sqrt(theta.ewens[2]+theta.s(segsites,seqlen,variance=T)[2])
              }
          }
      }
    retval
  }

