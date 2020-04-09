\name{landscape.new.epoch}
\alias{landscape.new.epoch}
\title{Create an Epoch}

\description{
  Create an epoch for a Rmetasim landscape object}
\usage{
  ## must be called AFTER integer, switch, and float params have
  ##   been created and after the demography has been created
  ## S, R, and M matricies must be square matricies of size X by X
  ##   where X = rland$intparam$stages*rland$intparam$habitats

rland <- landscape.new.epoch(rland,S=NULL,R=NULL,M=NULL,epochprob=1,
                               startgen=0,extinct=NULL,carry=NULL,localprob=NULL,
                               pollen.kernels = NULL, seed.kernels = NULL, 
                               leftx = NULL, rightx = NULL,
                               boty = NULL, topy = NULL, maxland = c(0, 0, 10000, 10000))
}

\arguments{
  
  \item{rland}{partially created landscape object, required}

  \item{S}{(default=NULL) Survivablity matrix for epoch, NULL gives no
    movement between subpopulations (0 matrix)}
  
  \item{R}{(default=NULL) female Reproduction matrix for epoch, NULL
    gives no dispersal between subpopulations other than that determined
    by dispersal kernels(0 matrix)}
  
  \item{M}{(default=NULL) Male reporduction matrix for epoch, NULL gives
    no sperm or pollen movement between subpopulations other then that
    determined by dispersal kernels (0 matrix)}
  
  \item{epochprob}{(default=1) probability of choosing this epoch
    randomly if randepoch==1}
  
  \item{startgen}{(default=0) generation in which this epoch starts}
  
  \item{extinct}{(default=NULL) vector of extinction probabilities per
    generation for each subpopulation, must be rland$intparam$habitats
    in length, passing NULL gives a 0\% probability of extinction to
    each subpopulation}
  
  \item{carry}{(default=NULL) vector of carrying capacities for each
    subpopulation, must be rland$intparam$habitats in length, passing NULL
    gives a 1000 individual carrying capacity to each subpopulation}
  
  \item{localprob}{(default=NULL) vector of probabilites for choosing
    local demographies, must be length(rland$demography$localdem) in
    length, passing NULL gives each demography an equalprobability}

  \item{leftx}{vector of the left x-coordinates of habitat patches.  If
  NULL (along with the other coordinates below), then the landscape is
  dividided into a set of contiguous strips}
  
  \item{rightx}{vector of the right x-coordinates of habitat patches}
  
  \item{topy}{vector of the bottom y-coordinates of habitat patches}

  \item{boty}{vector of the top y-coordinates of habitat patches}

  \item{pollen.kernels}{A matrix that describes pollen kernels for each
    stage in the landscape.  The rows correspond to stages in the
    landscape (there should be habitat * stages rows).  The columns
    correspond to the characteristics of the dispersal kernels: Column 1
    is the type of kernel (1=exponential, 2=Weibull, 3=mixture between
    Weibull and Gaussian) because the first is a special case of the
    second and the second is a special case of the third type of
    dispersal kernel using "3" in this column is pretty much always
    reasonable.  Columns 2 and 3 are scale and shape parameters of the
    Weibull (if shape=1, this reduces to an exponential kernel).
    Columns 4 and 5 correspond to the mean and sd of the Gaussian
    portion of the kernel.  Column 6 is the mixture parameter, ranging
    from 0-1 that gives the relative weights of the Weibull versus
    Gaussian portions of the kernel.  If the pollen kernel is constant
    across demographic stages, it is much easier to specify it when
    calling \code{landscape.new.floatparam()}}

    \item{seed.kernels}{A similar matrix to \code{pollen.kernels}}
}

  \examples{
  
  exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
  exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
  exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
  
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.intparam(exampleland, s=2, h=2)
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland <- landscape.new.switchparam(exampleland)
  exampleland <- landscape.new.local.demo(exampleland,exampleS,exampleR,exampleM)

  ## nonsense matricies
  exampleS <- matrix(c(rep(0,4),
                rep(1,4),
                rep(0,4),
                rep(1,4)), nrow = 4)
  exampleR <- matrix(c(rep(0.5,4),
                rep(0,4),
                rep(0.5,4),
                rep(0,4)), nrow = 4)
  exampleM <- matrix(c(rep(0,4),
                rep(.25,4),
                rep(0,4),
                rep(0,4)), nrow = 4)

  ## defaults
  exampleland<- landscape.new.epoch(exampleland,exampleS,exampleR,exampleM)

  exampleland$demography$epochs[[1]]

  rm(exampleS)
  rm(exampleR)
  rm(exampleM)
  rm(exampleland)

}
\keyword{misc}
