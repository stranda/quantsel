\name{landscape.genepop.output}
\alias{landscape.genepop.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ writes genepop files based on the genotypes of individuals in a landscape }
\description{
creates a genepop file.  Occupied habitats correspond to genepop populations.  
}
\usage{
landscape.genepop.output(l, fn = "genepop.out", title = "rmetasim landscape output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{landscape object}
  \item{fn}{output file name}
  \item{title}{title of output in genepop file}
}
\details{Only exports diploid data.  You might want to run \code{landscape.sample()} on the landscape before exporting it}

\seealso{ landscape.write.foreign, landscape.write.fdist}
\examples{
#	l <- landscape.new.example()	
#	l <- landscape.simulate(l,5)
#	landscape.genepop.output(landscape.sample(l,3,24))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}

