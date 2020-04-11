#
# formats a landscape locus to calculate amova
#




#' uses functions in ade4 to calcuate phi-st for a particular locus
#' 
#' Runs an amova on a locus.  Does not include information about sequence
#' similarity or ssr size in analysis.
#' 
#' Should be the same as Weir and Cockerham's theta
#' 
#' @param l locus number
#' @param rland landscape object
#' @return list of amova results for a locus
#' @seealso landscape.amova, landscape.amova.pairwise
#' @keywords misc
#' @export landscape.amova.locus
landscape.amova.locus <- function(l=1,rland)
{
  loc <- landscape.locus(l,rland)
  if (landscape.ploidy(rland)[l]>1)
    {
      sortloc <- t(apply(loc[,(landscape.democol()+1):(landscape.democol()+2)],1,sort))
      df <- data.frame(matrix(table(paste(sortloc[,1],sortloc[,2]),landscape.populations(rland)),ncol=length(unique(landscape.populations(rland))),byrow=F))
      names(df) <- colnames(table(paste(sortloc[,1],sortloc[,2]),landscape.populations(rland)))
      rownames(df) <- rownames(table(paste(sortloc[,1],sortloc[,2]),landscape.populations(rland)))
    }
  else
  {
    df <- data.frame(matrix(table(paste(loc[,landscape.democol()+1]),landscape.populations(rland)),
                            ncol=length(unique(landscape.populations(rland))),byrow=F))
    names(df) <- colnames(table(paste(loc[,landscape.democol()+1]),landscape.populations(rland)))
    rownames(df) <-   rownames(as.matrix(table(paste(loc[,landscape.democol()+1]),landscape.populations(rland))))
  }
  ade4::amova(samples=df)
}






#' calculates pairwise phi_ST for a landscape
#' 
#' pairwise phi_ST calculator.  Kind of slow. use \code{landscape.sample} to
#' reduce the size of the calculation.
#' 
#' 
#' @param rland landscape object
#' @seealso landscape.amova, landscape.amova.locus
#' @export landscape.amova.pairwise
landscape.amova.pairwise <- function(rland)
  {
    retmat <- matrix(NA,nrow=rland$intparam$habitat,ncol=rland$intparam$habitat)
    for (i in 1:rland$intparam$habitat)
      for (j in 1:i-1)
        {
          if (length(landscape.populations(rland))>0)
            {
              if ((max(landscape.populations(rland)==i)>0)&(max(landscape.populations(rland)==j)>0))
                {
                  compland <- rland
                  compland$individuals <- compland$individuals[landscape.populations(compland) %in% c(i,j),]
                  tmpvec <- rep(NA,length(rland$loci))
                  n <- 1
                  for (l in 1:length(rland$loci))
                    {
                      tmpvec[l] <- landscape.amova.locus(l,compland)$statphi
                      if (!is.na(tmpvec[l]))
                        n <- n+1;
                    }
             
                  retmat[i,j] <- sum(unlist(tmpvec),na.rm=T)/n
                }
            }
        }
    retmat
  }






#' calcuates phi-st for every locus in the landscape
#' 
#' calcuates phi_ST for every locus in the landscape
#' 
#' 
#' @param rland landscape object
#' @param np max number of pops to include
#' @param ns max number of samples to collect
#' @return vector of length equal to the number of loci
#' @seealso
#' \code{\link{landscape.amova.locus}},\code{\link{landscape.amova.pairwise}}
#' @keywords misc
#' @export landscape.amova
landscape.amova <- function(rland,np=24,ns=24) #np is the sample size per population.  Only pops with actual > np
  {
    unlist(sapply(1:length(rland$loci),function(x,l){landscape.amova.locus(x,l)$statphi
                                            },l=landscape.sample(rland,np=np,ns=ns)))
  }

