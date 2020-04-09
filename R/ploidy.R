"landscape.ploidy" <-
function(Rland)
  {
    ploidy<-c();
    for (i in 1:Rland$intparam$locusnum)
      {
        ploidy<-c(ploidy,Rland$loci[[i]]$ploidy);
      }
    ploidy
  }

