\name{landscape.generate.locations}
\alias{landscape.generate.locations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ generates randomly located and sized, non-overlapping habitat patches on the landscape}
\description{
	makes it easier to define landscapes, if you are just interested in a random landscape with certain qualities  
}
\usage{
landscape.generate.locations(npop = 10, xrange = c(0, 15000), yrange = xrange, sizexkernel = c(300, 80), sizeykernel = sizexkernel, boundaries = NULL)
}
\arguments{
  \item{npop}{number of habitats to define}
  \item{xrange}{ the low and high bounds for a landscape's x coordinates (stay positive) }
  \item{yrange}{ the low and high bounds for a landscape's y coordinates (stay positive) }
  \item{sizexkernel}{mean and sd for a normal distribution defining the x-dimension of habitats }
  \item{sizeykernel}{mean and sd for a normal distribution defining the y-dimension of habitats }
  \item{boundaries}{another way to define boundaries for landscape.}
}
\details{

Chooses x and y coordinates for rectangles from a uniform distribution
with limits given by the *range parameters.  The x and y sizes of the habitats
are governed by pulls from normal distributions

}
\value{
a matrix with columns 1=leftx x coordinates, 2=bottom y coordinates, 3=right x coordinates, 4=top y coordinates. The row number is equal to the number of habitats (\code{npop})
}

\seealso{ landscape.new.epoch()}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

