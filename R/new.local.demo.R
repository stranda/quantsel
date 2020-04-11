#' Create a Local Demography
#' 
#' Create a local demography for an Rmetasim Landscape object
#' 
#' The local demography objects encapsulate demography within a particular
#' region.  Multiple such objects can be defined to account for different
#' demographies across space.  The flag, k, can indicate whether the matrices
#' represent demography at zero population growth and at carrying capacity, if
#' density-dependence is modeled
#' 
#' @param rland partially created landscape object, required
#' @param S Survivablity matrix for demograpy, required
#' @param R female Reproduction matrix for demography, required
#' @param M Male reporduction matrix for demography, required
#' @param k flag for type of matrix, 0=demgraphy at zero population density,
#' 1=demography at carrying capacity
#' @keywords misc
#' @examples
#' 
#'   
#'   exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
#'   exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
#'   exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
#'   
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.intparam(exampleland, s=2)
#'   exampleland <- landscape.new.floatparam(exampleland)
#'   exampleland <- landscape.new.switchparam(exampleland)
#'   exampleland <- landscape.new.local.demo(exampleland,exampleS,exampleR,exampleM)
#' 
#'   exampleland$demography$localdem
#' 
#'   rm(exampleS)
#'   rm(exampleR)
#'   rm(exampleM)
#'   rm(exampleland)
#' 
#' @export landscape.new.local.demo
"landscape.new.local.demo" <-
function(rland,S,R,M)
{
  if (is.null(rland$demography$localdem))
    {
      rland$demography$localdem <- list(NULL)
      demonum <- 1
    }
  else
    {
      demonum <- length(rland$demography$localdem) + 1
      rland$demography$localdem[[demonum]] <- list(LocalS=NULL,LocalR=NULL,LocalM=NULL)
    }

  if (is.nsquare(S,rland$intparam$stages) && 
      is.nsquare(R,rland$intparam$stages) && 
      is.nsquare(M,rland$intparam$stages))
    {
      rland$demography$localdem[[demonum]]$LocalS <- S    
      rland$demography$localdem[[demonum]]$LocalR <- R    
      rland$demography$localdem[[demonum]]$LocalM <- M
      rland$intparam$numdemos <- length(rland$demography$localdem) 
    }
  else
    {
      stop("Matricies do not conform to stages set in intparam!")
    }
  rland
}

