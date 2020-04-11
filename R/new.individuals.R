#' Fill a landscape with individuals
#' 
#' Create a set of individuals for a Rmetasim landscape object.
#' 
#' 
#' @param rland nearly complete landscape object, required
#' @param PopulationSizes vector of integers denoting how many individuals are
#' in which stage and in which subpopulation, vector is ordered as: (pop1
#' stage1, pop1 stage2, ..., pop2 stage1, pop2stage2, ....), must be of length
#' rland$intparam$habitats * rland$intparam$stages
#' @keywords misc
#' @examples
#' 
#'   
#'   exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
#'   exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
#'   exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
#'   
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.intparam(exampleland, s=2, h=2)
#'   exampleland <- landscape.new.floatparam(exampleland)
#'   exampleland <- landscape.new.switchparam(exampleland)
#'   exampleland <- landscape.new.local.demo(exampleland,exampleS,exampleR,exampleM)
#' 
#'   ## nonsense matricies
#'   exampleS <- matrix(c(rep(0,4),
#'                 rep(1,4),
#'                 rep(0,4),
#'                 rep(1,4)), nrow = 4)
#'   exampleR <- matrix(c(rep(0.5,4),
#'                 rep(0,4),
#'                 rep(0.5,4),
#'                 rep(0,4)), nrow = 4)
#'   exampleM <- matrix(c(rep(0,4),
#'                 rep(.25,4),
#'                 rep(0,4),
#'                 rep(0,4)), nrow = 4)
#' 
#'   exampleland<- landscape.new.epoch(exampleland,exampleS,exampleR,exampleM)
#'   exampleland <- landscape.new.locus(exampleland,type=2,ploidy=2,mutationrate=.001,numalleles=5,allelesize=100)
#'   exampleland <- landscape.new.locus(exampleland,type=1,ploidy=1,mutationrate=.001,numalleles=3)
#'   exampleland <- landscape.new.locus(exampleland,type=0,ploidy=2,mutationrate=.004,numalleles=4)
#' 
#'   exampleland <- landscape.new.individuals(exampleland,
#'                  c(5,20,7,15))
#' 
#'   exampleland$individuals
#' 
#'   rm(exampleS)
#'   rm(exampleR)
#'   rm(exampleM)
#'   rm(exampleland)
#' 
#' @export landscape.new.individuals
"landscape.new.individuals" <-
function(rland, PopulationSizes)
  {
    if (!is.null(rland$intparam)&(rland$intparam$habitats*rland$intparam$stages==length(PopulationSizes)))
      rland <- .Call("populate_Rland",rland,PopulationSizes,PACKAGE="quantsel")
    else
      warning("either the Population Size vector does not match other landscape parameters or those parameters have not been declared")
    rland
  }

landscape.add.individuals <- function(rland, PopulationSizes)
  {
    if (!is.null(rland$intparam)&(rland$intparam$habitats*rland$intparam$stages==length(PopulationSizes)))
      {
        rland2 <- .Call("populate_Rland",rland,PopulationSizes,PACKAGE="quantsel")
        rland2$individuals[,3] <- rep(rland$intparam$currentgen,dim(rland2$individuals)[1])
      }
    else
      stop("either the Population Size vector does not match other landscape parameters or those parameters have not been declared")
    rland$individuals <- rbind(rland2$individuals,rland$individuals)
    rland$individuals <- rland$individuals[order(rland$individuals[,1]),]
    rland
  }

