

#' Create a Skeletal Landscape
#' 
#' Create a skeletal Rmetasim landscape ready to be configured
#' 
#' 
#' @aliases landscape.new.empty
#' @keywords misc
#' @export
#' @examples
#' 
#'   ## Only usage
#'   landscape.new.empty()


"landscape.new.empty" <-
function()
{
  tmpdemo <- list(localdem=NULL,epochs=NULL)
  rland <- list(intparam=NULL,switchparam=NULL,floatparam=NULL,demography=tmpdemo,loci=NULL,expression=NULL,gpmap=NULL,individuals=NULL)
  rland
}

