"landscape.locusvec" <-
function(Rland)
  {
    p<-landscape.ploidy(Rland);
    lv<-c();
    for (i in  1:Rland$intparam$locusnum)
      {
        lv<-c(lv,rep(i,p[i]));
      }
    lv
  }

