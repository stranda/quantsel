#
# export the data in fdist format
#
landscape.write.fdist <- function(l,popnum=20,samp=24)
  {
    if (is.landscape(l)) #input error check
      {
        fn="infile"
        
        cat("0\n",file=fn)
        l <- landscape.sample(l,np=popnum,ns=samp)
        cat(paste(length(unique(landscape.populations(l))),"\n"),file=fn,append=T)
        cat(paste(length(l$loci),"\n"),file=fn,append=T)

        for (i in 1:length(l$loci))
          {
            if (landscape.ploidy(l)[i]>1)
              {
                genotypes <- landscape.locus(i,l)[,-1:-landscape.democol()]
                alleles <- rbind(cbind(landscape.populations(l),genotypes[,1]),
                                 cbind(landscape.populations(l),genotypes[,1]))
                atbl <- table(alleles[,1],alleles[,2])
                cat(paste(dim(atbl)[2],"\n"),file=fn,append=T)
                tmp <- apply(atbl,1,function(x){cat(paste(x," ",sep=''),file=fn,append=T);cat("\n",file=fn,append=T)})
              }
          }
      }
  }

