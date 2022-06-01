#' Calculate allele frequencies at each locus in each population
#' 
#' Calculate allele frequencies
#' 
#' 
#' @param rland the Rmetasim landscape object
#' @param tbl.out return as a (three-dimensional) table if TRUE.  If FALSE,
#' return as a dataframe with categorical variables denoting the locus,
#' population and allele.
#' @return Depends on the value of tbl.out.  See above.
#' @seealso landscape.obs.het, landscape.exp.het, landscape.Fwright,
#' landscape.Fst
#' @keywords misc
#' @examples
#' 
#' #  exampleland <- landscape.new.example()
#' #  exampleland <- landscape.simulate(exampleland, 4)
#' #  landscape.allelefreq(exampleland,tbl.out=TRUE)
#' #  landscape.allelefreq(exampleland,tbl.out=FALSE)
#' #  rm(exampleland)
#' 
#' @export landscape.allelefreq
"landscape.allelefreq" <-
function(Rland,tbl.out=FALSE)
  {
    counts <- landscape.allelecount(Rland)
    popsizes <- data.frame(table(landscape.populations(Rland)))
    names(popsizes) <- c("pop","N")
    counts <- merge(counts,popsizes,all.x=T)
    ploidymult <- landscape.ploidy(Rland)[counts$loc]
    counts$Freq <- counts$Freq/(ploidymult*counts$N)
    if (tbl.out==T)
      {
        xtabs(Freq~pop+alleles+loc,counts)
      } else {
      counts[,1:4]
    }
  }

