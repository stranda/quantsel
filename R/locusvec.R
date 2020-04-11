#' return a vector with the locus ids for each column in the individuals
#' component of a landscape
#' 
#' return a vector with the locus ids for each column in the individuals
#' component of a landscape
#' 
#' 
#' @param Rland the Rmetasim landscape object
#' @return vector
#' @seealso landscape.populations
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   exampleland <- landscape.simulate(exampleland, 4)
#'   landscape.locusvec(exampleland)
#'   rm(exampleland)
#' 
#' @export landscape.locusvec
"landscape.locusvec" <-
function(Rland)
  {
    p<-landscape.ploidy(Rland);
    lv<-c();
    for (i in  1:Rland$intparam$locusnum)
      {
        lv<-c(lv,rep(i,p[i]));
      }
    lv
  }

