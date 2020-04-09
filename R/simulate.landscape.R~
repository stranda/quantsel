landscape.coerce <- function(rland)
  {
    rland$intparam <- lapply(rland$intparam,as.integer)
    rland$switchparam <- lapply(rland$switchparam,as.integer)
    rland$loci <- lapply(rland$loci,function(x)
                         {
                           x$type <- as.integer(x$type)
                           x$ploidy <- as.integer(x$ploidy)
                           x$trans <- as.integer(x$trans)
                           x$alleles <- lapply(x$alleles,function(y,typ)
                                               {
                                                 y$aindex <- as.integer(y$aindex)
                                                 y$birth <- as.integer(y$birth)
                                                 if (typ!=253)
                                                   {
                                                     y$state <- as.integer(y$state)
                                                   }
                                                 y
                                               },
                                               typ=x$type)
                           x
                         })
    rland$individuals <- matrix(as.integer(rland$individuals),nrow=dim(rland$individuals)[1])
    rland
  }

#this function implements carry capacity only in stage 0 of each habitat
landscape.simulate.stg0 <- function(Rland, numit, seed=-1, compress=FALSE)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("iterate_landscape_stg0",as.integer(numit),Rland,as.integer(compress),PACKAGE = "kernelPop2")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }


landscape.simulate <- function(Rland, numit, seed=-1, compress=FALSE)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("iterate_landscape",as.integer(numit),Rland,as.integer(compress),PACKAGE = "kernelPop2")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

