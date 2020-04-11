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
        .Call("iterate_landscape_stg0",as.integer(numit),Rland,as.integer(compress),PACKAGE = "quantsel")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }






#' Run a simulation for a single landscape through time
#' 
#' Simulate a Rmetasim landscape for a number of generations.
#' 
#' 
#' @param rland the Rmetasim landscape object
#' @param numit the number of generations/iterations to simulate, note that
#' landscapes will not run past the rland$intparam$totalgens value
#' @param seed The default value of seed uses the seed set in the calling
#' environment.  Any other value for seed uses 'set.seed()' to reset the random
#' number generator.  landscape.simulate uses the RNG selected by the calling
#' environment.
#' @param compress If true, landscape.simulate executes a survival and carrying
#' capacity step before returning.  In demographies with high reproductive
#' potential, this can significantly reduce the size of R objects returned
#' @param adj.lambda Tries to apply a correction to population growth that
#' makes the observed growth rate more closely approximate that predicted from
#' standard analysis eigensystem of the sum of the survival and reproduction
#' Lefkovitch matrices
#' @keywords misc
#' @examples
#' 
#' #  exampleland <- landscape.new.example()
#' #  exampleland <- landscape.simulate(exampleland, 4)
#' #  exampleland
#' #  rm(exampleland)
#' 
#' @export landscape.simulate
landscape.simulate <- function(Rland, numit, seed=-1, compress=FALSE)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("iterate_landscape",as.integer(numit),Rland,as.integer(compress),PACKAGE = "quantsel")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

