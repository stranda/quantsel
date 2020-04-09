#
#
# A quick R function to take a landscape and output a genepop file
# AES 1/26/06
#
# This function will just ignore haploid loci and only output diploids.  It also uses the three digit allele designations
# It takes all individuals and puts them into the output file
#
# IMPORTANT to avoid the missing data indicator in GENEPOP, 1 is added to all allele numbers
#
# l is a landscape and fn is the name of the output file
#
landscape.genepop.output <- function(l,fn="genepop.out",title="rmetasim landscape output")
  {
    if (is.landscape(l)) #input error check
      {
        cat(paste(date(),": ",title,sep=''),file=fn)
        cat("\n",file=fn,append=T)
        #first output the locus names
        for (i in which(landscape.ploidy(l)==2)) #only use the diploids
          {
            cat(paste("locus-",i," ",sep=''),file=fn,append=T)
            cat("\n",file=fn,append=T)
          }
        #now outputting the genotypic data is a bit more complex
        #first make a big matrix of the genotypes is GENEPOP format
        #the rows still correspond to the rows in l$individuals
        printmat <- matrix("",ncol=(1+length(which(landscape.ploidy(l)==2))+1),
                           nrow=length(landscape.populations(l)))
        printmat[,1] <- sprintf("Class-%d Gen-%d,",l$individuals[,1],l$individuals[,3])
        colcount <- 2
        for (i in which(landscape.ploidy(l)==2))
          {
            if (l$loci[[i]]$type!=253)
              {
                printmat[,colcount] <- sprintf("%03d%03d ",landscape.states(i,l)[,(landscape.democol()+1)]+1,
                                               landscape.states(i,l)[,(landscape.democol()+2)]+1)
              }
            else
              {
                printmat[,colcount] <- sprintf("%03d%03d ",landscape.locus(i,l)[,(landscape.democol()+1)]+1,
                                               landscape.locus(i,l)[,(landscape.democol()+2)]+1)
              }
            colcount <- colcount + 1
          }
        printmat[,colcount] <- rep("\n",length(landscape.populations(l)))
        #now output one population at a time
        for (i in unique(landscape.populations(l)))
          {
            cat("POP\n",file=fn,append=T)
#            apply(printmat[populations(l)==i,],1,function(x,fn) {cat (x, file=fn,append=T)},fn=fn)
            cat(t(printmat[landscape.populations(l)==i,]),file=fn,append=T)
          }
      } else { #see, it really wasn't a landscape
        print ("make sure that the landscape is in the right form")
      }
  }
