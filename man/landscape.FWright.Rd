\name{landscape.FWright}
\alias{landscape.FWright}
\title{Calculates Wright's Fixation index at each locus in each population}

\description{
  Calculate F
  }
\usage{
  Fmat <- landscape.FWright(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\value{
  A matrix of population x locus measures of Wright's F
}
\seealso{
obs.het.landscape, exp.het.landscape, allelefreq.landscape, Fst.landscape
}
\examples{
  exampleland <- landscape.new.example()
  landscape.FWright(exampleland)
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.FWright(exampleland)
  rm(exampleland)
}
\keyword{misc}
