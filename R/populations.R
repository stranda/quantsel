#' return a vector of population IDs from a landscape
#' 
#' return a vector of population IDs from a landscape
#' 
#' Returns a vector of length \code{dim(rland$individuals)[1]} where rland is a
#' landscape object.  The vector classifies individuals into populations (or
#' habitats)
#' 
#' @param Rland the Rmetasim landscape object
#' @return a vector
#' @seealso landscape.locus, landscape.ploidy
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   exampleland <- landscape.simulate(exampleland, 4)
#'   plot(table(landscape.populations(exampleland)),main="Distribution of population size in landscape")
#'   rm(exampleland)
#' 
#' @export landscape.populations
"landscape.populations" <-
function(Rland)
  {
    if (is.landscape(Rland,verb=F))
      {
        (Rland$individuals[,1]%/%Rland$intparam$stages)+1
      }
  }

