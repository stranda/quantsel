
"landscape.locus" <-
function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          Rland$individuals[,c(rep(TRUE,landscape.democol()),landscape.locusvec(Rland)==lnum)]
        }
  }


#returns a individual x ploidy matrix of states
landscape.states.old <- function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          lmat <- as.data.frame(landscape.locus(lnum,Rland))
          st <- landscape.states(lnum,Rland)
          stmat <- matrix(NA,ncol=2,nrow=max(st$aindex)+1)
          stmat[st$aindex+1,1] <- st$aindex
          stmat[st$aindex+1,2] <- st$state
          
          lmat[,landscape.democol()+1] <- stmat[lmat[,landscape.democol()+1]+1,2]

          if (landscape.ploidy(Rland)[lnum]==2)
            {
              lmat[,landscape.democol()+2] <- stmat[lmat[,landscape.democol()+2]+1,2]
            }
          lmat
        }
  }

#
#takes a locus and returns the states and their indices
#
landscape.locus.states<-function(lnum=1,Rland)
  {
    if (is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          ain<-c();
          sta<-c();
          locin <- landscape.locus(lnum,Rland)[,c(-1:-landscape.democol())]
#          print(locin)
          ainds <- unique(c(locin))
#          print(ainds)
          for (i in 1:length(Rland$loci[[lnum]]$alleles))
            {
              if (Rland$loci[[lnum]]$alleles[[i]]$aindex %in% ainds)
                {
                  ain<-c(ain,Rland$loci[[lnum]]$alleles[[i]]$aindex);
                  sta<-c(sta,Rland$loci[[lnum]]$alleles[[i]]$state);
                }
            }
          list(aindex=ain,state=sta);
          
        }
  }
