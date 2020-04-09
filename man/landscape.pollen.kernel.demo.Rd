\name{landscape.pollen.kernel.demo}
\alias{landscape.pollen.kernel.demo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ takes parameters that define a pollen kernel and plots the results of simulating that kernel}
\description{
Produces a plot of a landscape.  Adds a plot of the pollination distances in the landscape.
}
\usage{
landscape.pollen.kernel.demo(s = 0, pmn1 = 30, pshp1 = 1.01, pmn2 = 600, pshp2 = 200, pmix = 0.75, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{ selfing rate to define}
  \item{pmn1}{ mean of component distribution 1 }
  \item{pshp1}{ shape of component distribution 1 }
  \item{pmn2}{mean of component distribution 2 }
  \item{pshp2}{ shape of component distribution 2  }
  \item{pmix}{ mixing parameter [0,1], 1=> all component 1}
  \item{plot}{ boolean }
}
\value{
landscape object 
}

\examples{
l<-landscape.pollen.kernel.demo()
}

\keyword{ misc }% __ONLY ONE__ keyword per line
