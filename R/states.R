"landscape.states" <-
function(lnum=1,Rland)
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
          sl <- list(aindex=ain,state=sta);
          if (landscape.ploidy(Rland)[lnum]>1)
            {
              st <- cbind(sl$state[match(locin[,1],sl$aindex)],sl$state[match(locin[,2],sl$aindex)])
            }
          else
            {
              st <- matrix(sl$state[match(c(locin),sl$aindex)],ncol=1)
            }
          cbind(landscape.locus(lnum,Rland)[,1:landscape.democol()],st)
        }
    
    
  }

