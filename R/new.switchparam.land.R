#' Create a set of boolean parameters
#' 
#' Create a set of boolean (1 or 0) parameters for a Rmetasim landscape.
#' 
#' 
#' @param rland skeletion of landscape object, required
#' @param re randepoch (default=0), 1=randomly pick a new epoch (from the
#' epochs listed in the landscape) after an epoch completes, 0=epochs are
#' chosen in order
#' @param rd randdemo (default=0), 1=randomly choose a demography (from the
#' demographies listed in the landscape) for each subpopulation, 0=demographies
#' are assigned in order
#' @param mp multp (default=1), 1=multiple paternity,0=entire families from a
#' single mating
#' @param dd density dependence.  If dd=1, then two of each local demography
#' matrix must be defined, the first set using new.local.demo with k=0 and
#' representing demography at low density and again with k=1 for demography at
#' high population density.
#' @keywords misc
#' @examples
#' 
#'   ## Defaults
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.switchparam(exampleland)
#'   exampleland$switchparam
#' 
#'   ## Random epochs, random demographies, and no multiple paternity
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.switchparam(exampleland,re=1,rd=1,mp=0)
#'   exampleland$switchparam
#' 
#'   rm(exampleland)
#' 
#' @export landscape.new.switchparam
"landscape.new.switchparam" <-
function(rland, re=0,rd=0,mp=1,sc=1)
{
  rland$switchparam <- list(re,rd,mp)
  names(rland$switchparam) <- c( "randepoch","randdemo","multp")
  rland
}

