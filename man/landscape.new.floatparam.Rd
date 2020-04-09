\name{landscape.new.floatparam}
\alias{landscape.new.floatparam}
\title{Create a set of floating point parameters}

\description{
  Create a set of floating point parameters for a Rmetasim landscape.}
\usage{
  ## must be called AFTER landscape.new.empty()
  rland <- landscape.new.floatparam(rland,s=0, seedscale = c(10,10),
  seedshape = c(10,10), seedmix=1, pollenscale=c(2,10),
  pollenshape=c(2,6), pollenmix=1,asp=1)}

\arguments{
  \item{rland}{skeletion of landscape object, required}

  \item{s}{selfing (default=0), the selfing rate of the species}


  \item{seedscale}{two element vector that provides the scale of the
  exponential/Weibull and mean of the Gaussian components of the seed
  dispersal kernel.}
  
  \item{seedshape}{two element vector that provides the shape of the
  exponential/Weibull and sd of the Gaussian components of the seed
  dispersal kernel, respectively.  If the first element of the vector is
  set to 1, the Weibull portion of the kernel becomes an exponential.}

  \item{seedmix}{the mixture parameter (ranges from 0-1) that weights
  the Weibull/exponential versus Gaussian portions of the dispersal
  kernel.  A value of 1 selects dispersal distances from only the
  Weibull/exponential portion.  A value of 0 selects distances only from
  the Gaussian portion.}


  \item{pollenscale}{two element vector that provides the scale of the
  exponential/Weibull and mean of the Gaussian components of the pollen
  dispersal kernel.}
  
  \item{pollenshape}{two element vector that provides the shape of the
  exponential/Weibull and sd of the Gaussian components of the pollen
  dispersal kernel, respectively.  If the first element of the vector is
  set to 1, the Weibull portion of the kernel becomes an exponential.}

  \item{pollenmix}{the mixture parameter (ranges from 0-1) that weights
  the Weibull/exponential versus Gaussian portions of the dispersal
  kernel.  A value of 1 selects dispersal distances from only the
  Weibull/exponential portion.  A value of 0 selects distances only from
  the Gaussian portion.}

}
\examples{
  ## Defaults
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland$floatparam

  ## .5 selfing rate
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.floatparam(exampleland,s=0.5)
  exampleland$floatparam

  rm(exampleland)
}
\keyword{misc}
