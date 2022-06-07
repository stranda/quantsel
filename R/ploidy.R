#' return a vector with the ploidy of each locus
#' 
#' return a vector with the ploidy of each locus in the order they appear in
#' the landscape
#' 
#' 
#' @param Rland the Rmetasim landscape object
#' @return vector
#' @seealso landscape.populations
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   landscape.ploidy(exampleland)
#'   rm(exampleland)
#' 
#' @export landscape.ploidy
"landscape.ploidy" <-
function(Rland)
  {
    ploidy<-c();
    for (i in 1:Rland$intparam$locusnum)
      {
        ploidy<-c(ploidy,Rland$loci[[i]]$ploidy);
      }
    ploidy
  }

