\name{landscape.write.fdist}
\alias{landscape.write.fdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writes a landscapes allele frequencies to an FDIST file}
\description{
Calculates the allele frequencies in each occupied habitat up to the number specified in popnum.  Makes calculations based upon the number of individuals sampled stored in samp
}
\usage{
landscape.write.fdist(l, popnum = 20, samp = 24)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{ landscape object}
  \item{popnum}{ max number of populations to sample }
  \item{samp}{ max number of individuals to sample per population}
}
\value{
	NULL
}

\seealso{\code{\link{landscape.write.foreign}},\code{\link{landscape.genepop.output}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

