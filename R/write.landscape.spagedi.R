#
#
# outputs a landscape in SPAgEDi
# if you want to subset the landscape, you need to do so 
# before calling this program
# only writes out the diploids and basically ignores the population
# designations, just focuses on the spatial coordinates
#
#

landscape.write.spagedi <- function(l,dists=10,
                                    fn="spagedi.txt",
                                    title="kernelPop output")
                                        #l is a kernelpop landscape
  {
    if (!is.landscape(l,verb=F)) stop ("must use a valid landscape")
###title
    cat(paste("//",date(),": ",title,sep=''),file=fn)
    cat("\n",file=fn,append=T)
###first line
    cat(dim(l$ind)[1],file=fn,append=T)
    cat("\t0\t2\t",file=fn,append=T)
    cat(length(which(landscape.ploidy(l)==2)),file=fn,append=T)
    cat("\t3\t2",file=fn,append=T)
    cat("\n",file=fn,append=T)
###second line
    if (length(dists)>1)
      {
        cat(paste(dists,"\t",sep=""),file=fn,append=T)
      } else
    {
      cat(paste("-",dists,sep=""),file=fn,append=T)
    }
    cat("\n",file=fn,append=T)
###third line
    cat(paste("Ind","X","Y",sep="\t"),file=fn,append=T)
    for (i in 1:length(which(landscape.ploidy(l)==2)))
      {
        cat("\t",file=fn,append=T)
        cat(paste("locus",i,sep=""),file=fn,append=T)
      }
    cat("\n",file=fn,append=T)
    
####now output the genotypic data

    
    ##now outputting the genotypic data is a bit more complex
    ##first make a big matrix of the genotypes is GENEPOP format
    ##the rows still correspond to the rows in l$individuals
    printmat <- matrix("",ncol=(1+length(which(landscape.ploidy(l)==2))+3),
                       nrow=length(landscape.populations(l)))
    printmat[,1] <- sprintf("stg%d-Gen%d\t",l$individuals[,1],l$individuals[,3])
    printmat[,2] <- sprintf("%d\t",l$individuals[,4])
    printmat[,3] <- sprintf("%d\t",l$individuals[,5])
    colcount <- 4
    for (i in which(landscape.ploidy(l)==2))
      {
        if (l$loci[[i]]$type!=253)
          {
            printmat[,colcount] <-
              sprintf("%03d%03d\t",landscape.states(i,l)[,(landscape.democol()+1)]+1,
                      landscape.states(i,l)[,(landscape.democol()+2)]+1)
          }
        else
          {
            printmat[,colcount] <-
              sprintf("%03d%03dt",landscape.locus(i,l)[,(landscape.democol()+1)]+1,
                      landscape.locus(i,l)[,(landscape.democol()+2)]+1)
          }
        colcount <- colcount + 1
      }
    printmat[,colcount] <- rep("\n",length(landscape.populations(l)))
    cat(t(printmat),file=fn,append=T)
    cat("END",file=fn,append=T)
    cat("\n",file=fn,append=T)
    
  }
