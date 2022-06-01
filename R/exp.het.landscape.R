#' Calculate expected heterozygosity
#' 
#' Calculate expected heterozygosity from a landscape
#' 
#' Calculates the expected heterozygosity in each population:
#' \deqn{1-\Sigma_{i_k} p_i^2}{1 - sum(p^2)} where \eqn{p} is a vector of
#' allele frequencies for a locus in a population.
#' 
#' @param rland the Rmetasim landscape object
#' @return A matrix with num loci columns and num populations rows.  Each
#' element reflects the expected heterozygosity for that population x locus
#' combination
#' @seealso landscape.obs.het, Fst.landscape
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   exampleland <- landscape.simulate(exampleland, 4)
#'   exphet <- landscape.exp.het(exampleland)
#'   rm(exampleland)
#' 
#' @export landscape.exp.het
landscape.exp.het <- function (Rland) 
{
    pl <- landscape.ploidy(Rland)
    strt <-  9 + cumsum(pl) - pl + 1  
    stp <- strt + (pl - 1)
    rl <- matrix(0, nrow = Rland$intparam$habitats, ncol = length(Rland$loci))
    afrq <- landscape.allelefreq(Rland)
    homframe <- aggregate(afrq$Freq^2,by=list(pop=afrq$pop,loc=afrq$loc),sum)
    homframe$het <- 1-homframe$x
    homframe <- homframe[,-3]
    as.matrix(reshape(homframe,idvar="pop",timevar="loc",direction="wide")[,-1])
  }
