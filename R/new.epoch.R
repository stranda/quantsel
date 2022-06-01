#' Create an Epoch
#' 
#' Create an epoch for a Rmetasim landscape object
#' 
#' 
#' @param rland partially created landscape object, required
#' @param S (default=NULL) Survivablity matrix for epoch, NULL gives no
#' movement between subpopulations (0 matrix)
#' @param R (default=NULL) female Reproduction matrix for epoch, NULL gives no
#' dispersal between subpopulations other than that determined by dispersal
#' kernels(0 matrix)
#' @param M (default=NULL) Male reporduction matrix for epoch, NULL gives no
#' sperm or pollen movement between subpopulations other then that determined
#' by dispersal kernels (0 matrix)
#' @param epochprob (default=1) probability of choosing this epoch randomly if
#' randepoch==1
#' @param startgen (default=0) generation in which this epoch starts
#' @param extinct (default=NULL) vector of extinction probabilities per
#' generation for each subpopulation, must be rland$intparam$habitats in
#' length, passing NULL gives a 0\% probability of extinction to each
#' subpopulation
#' @param carry (default=NULL) vector of carrying capacities for each
#' subpopulation, must be rland$intparam$habitats in length, passing NULL gives
#' a 1000 individual carrying capacity to each subpopulation
#' @param localprob (default=NULL) vector of probabilites for choosing local
#' demographies, must be length(rland$demography$localdem) in length, passing
#' NULL gives each demography an equalprobability
#' @param leftx vector of the left x-coordinates of habitat patches.  If NULL
#' (along with the other coordinates below), then the landscape is dividided
#' into a set of contiguous strips
#' @param rightx vector of the right x-coordinates of habitat patches
#' @param topy vector of the bottom y-coordinates of habitat patches
#' @param boty vector of the top y-coordinates of habitat patches
#' @param pollen.kernels A matrix that describes pollen kernels for each stage
#' in the landscape.  The rows correspond to stages in the landscape (there
#' should be habitat * stages rows).  The columns correspond to the
#' characteristics of the dispersal kernels: Column 1 is the type of kernel
#' (1=exponential, 2=Weibull, 3=mixture between Weibull and Gaussian) because
#' the first is a special case of the second and the second is a special case
#' of the third type of dispersal kernel using "3" in this column is pretty
#' much always reasonable.  Columns 2 and 3 are scale and shape parameters of
#' the Weibull (if shape=1, this reduces to an exponential kernel).  Columns 4
#' and 5 correspond to the mean and sd of the Gaussian portion of the kernel.
#' Column 6 is the mixture parameter, ranging from 0-1 that gives the relative
#' weights of the Weibull versus Gaussian portions of the kernel.  If the
#' pollen kernel is constant across demographic stages, it is much easier to
#' specify it when calling \code{landscape.new.floatparam()}
#' @param seed.kernels A similar matrix to \code{pollen.kernels}
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
#'   ## defaults
#'   exampleland<- landscape.new.epoch(exampleland,exampleS,exampleR,exampleM)
#' 
#'   exampleland$demography$epochs[[1]]
#' 
#'   rm(exampleS)
#'   rm(exampleR)
#'   rm(exampleM)
#'   rm(exampleland)
#' 
#' 
#' @export landscape.new.epoch
"landscape.new.epoch" <- function(rland,S=NULL,R=NULL,M=NULL,epochprob=1,startgen=0,extinct=NULL,carry=NULL,localprob=NULL,
                                  pollen.kernels=NULL,
                                  seed.kernels=NULL,
                                  leftx=NULL, rightx=NULL, boty=NULL, topy=NULL,
                                  maxland=c(0,0,10000,10000))
{
  popw=0.1 * sum(maxland[c(1,3)])/2
  poph=0.1 * sum(maxland[c(2,4)])/2
  if (is.null(rland$demography$epochs))
    {
      rland$demography$epochs <- list(NULL)
      epochnum <- 1
    }
  else
    {
      epochnum <- length(rland$demography$epochs) + 1
      rland$demography$epochs[[epochnum]] <- list(RndChooseProb=NULL,StartGen=NULL,Extinct=NULL,
                                                  Carry=NULL,Localprob=NULL,S=NULL,R=NULL,M=NULL)
    }
  
  rland$demography$epochs[[epochnum]]$RndChooseProb <- epochprob
  rland$demography$epochs[[epochnum]]$StartGen <- startgen
  if (is.null(extinct))
    {
      extinct <- rep(0, rland$intparam$habitats)
    }  
  if (length(extinct) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$Extinct <- extinct
    }
  else
    {
      stop("Wrong size for extinct vector", dim((extinct)))
    }
  
  if (is.null(carry)) 
    {
      carry <- rep(1000, rland$intparam$habitats)
    }
  if (length(carry) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$Carry <- carry
    }
  else
    {
      stop("Wrong size for carry vector")
    }

  numdem <- length(rland$demography$localdem)
  if (is.null(localprob)) 
    {
      localprob <- rep(1/numdem, numdem)
    }
  if (length(localprob == numdem))
    {
      rland$demography$epochs[[epochnum]]$Localprob <- localprob
    }
  else
    {
      stop("Wrong size for localprob vector")
    }
  
  if (is.null(leftx)|is.null(rightx)|is.null(topy)|is.null(boty))
    {
      leftx <- rep(seq(maxland[1],maxland[3]-popw,length=ceiling(rland$intparam$habitats/2)),2)
      if (length(leftx)>rland$intparam$habitats)
        {
          leftx <- leftx[-length(leftx)]
        }
      rightx <- leftx+popw
      boty <- rep(seq(maxland[2],maxland[4]-poph,length=ceiling(rland$intparam$habitats/2)),2)
      if (length(leftx)>rland$intparam$habitats)
        {
          boty <- boty[-length(boty)]
        }
      topy <- boty+poph
    }  
  if (length(leftx) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$leftx  <- leftx
      rland$demography$epochs[[epochnum]]$rightx <- rightx
      rland$demography$epochs[[epochnum]]$boty   <- boty
      rland$demography$epochs[[epochnum]]$topy   <- topy
    }
  else
    {
      stop("Wrong size for coords vector", dim((extinct)))
    }

  matrixsize <- rland$intparam$habitats * rland$intparam$stages
  
  if (is.null(S))
    {
      S <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.null(R))
    {
      R <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.null(M))
    {
      M <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.nsquare(S,matrixsize) && 
      is.nsquare(R,matrixsize) && 
      is.nsquare(M,matrixsize))
    {
      rland$demography$epochs[[epochnum]]$S <- S
      rland$demography$epochs[[epochnum]]$R <- R
      rland$demography$epochs[[epochnum]]$M <- M
      rland$intparam$numepochs <- length(rland$demography$epochs)
    }
  else
    {
      stop("S, R, and M matricies are not the correct size")
    }

  if(is.null(seed.kernels))
    {# a kernel matrix has the same number of rows as a side of the landscape matrix and
      #six columns.  The columns correspond to 1 kernel, 2-3 scale and shape parameters of the kernel
      #4-5 scale and shape parameters for the 2nd distribution in a mixture and 6 mixing parameter (proportion of the first dist)
      seed.kernels <- cbind(rep(3,matrixsize),rep(rland$floatparam$seedmu,matrixsize),
                            rep(rland$floatparam$seedshape,matrixsize),rep(rland$floatparam$seedmu2,matrixsize),
                            rep(rland$floatparam$seedshape2,matrixsize),rep(rland$floatparam$seedmix,matrixsize))
    }
  if ((dim(seed.kernels)[2]!=6)||(dim(seed.kernels)[1]!=rland$intparam$habitats*rland$intparam$stages))
    {
      stop("Kernel matrix of incorrect dimension or landscape.intparam not run yet")
    }

  if(is.null(pollen.kernels))
    {# a kernel matrix has the same number of rows as a side of the landscape matrix and
      #six columns.  The columns correspond to 1 kernel, 2-3 scale and shape parameters of the kernel
      #4-5 scale and shape parameters for the 2nd distribution in a mixture and 6 mixing parameter (proportion of the first dist)
      pollen.kernels <- cbind(rep(3,matrixsize),rep(rland$floatparam$pollenmu,matrixsize),rep(rland$floatparam$pollenshape,matrixsize),rep(rland$floatparam$pollenmu2,matrixsize),rep(rland$floatparam$pollenshape2,matrixsize),rep(rland$floatparam$pollenmix,matrixsize))
    }
  
  if ((dim(pollen.kernels)[2]!=6)||(dim(pollen.kernels)[1]!=rland$intparam$habitats*rland$intparam$stages))
    {
      stop("Kernel matrix of incorrect dimension or landscape.intparam not run yet")
    }
  
  rland$demography$epochs[[epochnum]]$seedkern <- seed.kernels
  rland$demography$epochs[[epochnum]]$pollenkern <- pollen.kernels
  rland
}

