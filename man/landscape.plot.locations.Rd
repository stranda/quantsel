\name{landscape.plot.locations}
\alias{landscape.plot.locations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plots the locations of all habitats and individuals }
\description{

Plots the landscape.  Gives different habitats different colors,
though the colors cycle through rapidly when the number of habitats is
large.  Offspring always have the color of their \emph{mother's}
habitat.  This provides a quick way to assess inter-habitat movement
visually

}
\usage{
landscape.plot.locations(rland, adults = c(NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rland}{landscape}
  \item{adults}{vector of stages to select and plot, if NULL, plot all stages}
}

\value{
NULL
}

\keyword{misc}

