#
# phenotypes assuming additivity ans expression matrix control
#
# these phenotypes are just the 'multipliers' used in the C++ code
#
landscape.phenotypes <- function(rland)
  {
    rv <- matrix(0,nrow=dim(rland$individuals)[1],ncol=rland$intparam$nphen)
    
    for (p in 1:rland$intparam$nphen)
      {
        ac <- rep(0,dim(rland$individuals)[1])
        for (l in 1:length(rland$loci))
          {
            if (rland$loci[[l]]$type %in% (252))
              {
                st <- as.matrix(landscape.states(l,rland)[,-1:-9])
                ac <- ac + st[,1]*rland$expression$expmat[l,p]
                if (landscape.ploidy(rland)[l]>1)
                  ac <- ac+st[,2]*rland$expression$expmat[l,p]

              }
          }
        rv[,p] <- ac + rnorm(length(ac),mean=0,sd=mean(ac)*(1-rland$expression$hsq[p]))
        rv[rv[,p]<0,p] <- 0
      }
    rv
  }


landscape.phenotypes.c <- function(rland)
{

    rland$individuals=as.matrix(rland$individuals)
    
    matrix(.Call("phenotypes",rland,PACKAGE = "kernelPop2"),
           ncol=rland$intparam$nphen,byrow=T)

}

