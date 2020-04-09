
"landscape.allelefreq" <-
function(Rland,tbl.out=FALSE)
  {
    counts <- landscape.allelecount(Rland)
    popsizes <- data.frame(table(landscape.populations(Rland)))
    names(popsizes) <- c("pop","N")
    counts <- merge(counts,popsizes,all.x=T)
    ploidymult <- landscape.ploidy(Rland)[counts$loc]
    counts$Freq <- counts$Freq/(ploidymult*counts$N)
    if (tbl.out==T)
      {
        xtabs(Freq~pop+alleles+loc,counts)
      } else {
      counts[,1:4]
    }
  }

