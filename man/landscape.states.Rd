\name{landscape.states}
\alias{landscape.states}

\title{return a list object  containing actual allele states and their associated indices for a particular locus}
\description{return a list object  containing actual allele states and their associated indices for a particular locus}

\usage{
  landscape.states(lnum=1,Rland)
}

\arguments{
  \item{lnum}{the locus to return}
  \item{Rland}{the Rmetasim landscape object}
}


\value{
  list 
}

\details{
Returns a list with two elements, aindex containing the
allele indices and state containing the actual allele states.
}

\seealso{
  landscape.locus, landscape.locus.states
  }

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.states(1,exampleland)
  rm(exampleland)
}

\keyword{misc}
