#probably needs a good bit of error checking written in
#




#' Create a set of floating point parameters
#' 
#' Create a set of floating point parameters for a Rmetasim landscape.
#' 
#' 
#' @param rland skeletion of landscape object, required
#' @param s selfing (default=0), the selfing rate of the species
#' @param seedscale two element vector that provides the scale of the
#' exponential/Weibull and mean of the Gaussian components of the seed
#' dispersal kernel.
#' @param seedshape two element vector that provides the shape of the
#' exponential/Weibull and sd of the Gaussian components of the seed dispersal
#' kernel, respectively.  If the first element of the vector is set to 1, the
#' Weibull portion of the kernel becomes an exponential.
#' @param seedmix the mixture parameter (ranges from 0-1) that weights the
#' Weibull/exponential versus Gaussian portions of the dispersal kernel.  A
#' value of 1 selects dispersal distances from only the Weibull/exponential
#' portion.  A value of 0 selects distances only from the Gaussian portion.
#' @param pollenscale two element vector that provides the scale of the
#' exponential/Weibull and mean of the Gaussian components of the pollen
#' dispersal kernel.
#' @param pollenshape two element vector that provides the shape of the
#' exponential/Weibull and sd of the Gaussian components of the pollen
#' dispersal kernel, respectively.  If the first element of the vector is set
#' to 1, the Weibull portion of the kernel becomes an exponential.
#' @param pollenmix the mixture parameter (ranges from 0-1) that weights the
#' Weibull/exponential versus Gaussian portions of the dispersal kernel.  A
#' value of 1 selects dispersal distances from only the Weibull/exponential
#' portion.  A value of 0 selects distances only from the Gaussian portion.
#' @keywords misc
#' @examples
#' 
#'   ## Defaults
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.floatparam(exampleland)
#'   exampleland$floatparam
#' 
#'   ## .5 selfing rate
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.floatparam(exampleland,s=0.5)
#'   exampleland$floatparam
#' 
#'   rm(exampleland)
#' 
#' @export landscape.new.floatparam
landscape.new.floatparam <-
function(rland,s=0,seedscale=c(10,10),seedshape=c(10,10),seedmix=1,
         pollenscale=c(2,10),pollenshape=c(2,6),pollenmix=1,asp=1,
         mindens=1e-25)
{
  rland$floatparam <- list(selfing=s,seedmu=seedscale[1],seedshape=seedshape[1],
                           seedmu2=seedscale[2],seedshape2=seedshape[2],seedmix=seedmix[1],
                           pollenmu=pollenscale[1],pollenshape=pollenshape[1],
                           pollenmu2=pollenscale[2],pollenshape2=pollenshape[2],
                           pollenmix=pollenmix,aspect=asp,mindens=mindens)
  rland
}

