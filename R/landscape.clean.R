#' Function to resolve inconsistencies within a landscape
#' 
#' Converts a landscape to internal format and back.  This can resolve
#' inconsistencies in a 'hand-built' landscape
#' 
#' 
#' @param rland the Rmetasim landscape object
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   exampleland <- landscape.simulate(exampleland, 4)
#'   exampleland.clean <- landscape.clean(exampleland)
#'   rm(exampleland)
#' 
#' @export landscape.clean
landscape.clean <- function(rland)
  {
    if (is.landscape(rland))
      .Call("clean_landscape",rland,PACKAGE = "quantsel")
  }
