"indxfreq" <-
function(lnum=1,Rland)
  {
    lv<-landscape.locus(lnum,Rland)[,c(rep(FALSE,landscape.democol()),rep(TRUE,(ncol(landscape.locus(lnum,Rland))-landscape.democol())))];
    if (landscape.ploidy(Rland)[lnum]==1)
      {
        table(landscape.populations(Rland),lv)
      }
    else
      {
        lv2<-c(lv[,1],lv[,2]);
        table(rep(landscape.populations(Rland),2),lv2)
      }
    
  }

