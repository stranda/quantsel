#' return a list object containing actual allele states and their associated
#' indices for a particular locus
#' 
#' return a list object containing actual allele states and their associated
#' indices for a particular locus
#' 
#' Returns a list with two elements, aindex containing the allele indices and
#' state containing the actual allele states.
#' 
#' @param lnum the locus to return
#' @param Rland the Rmetasim landscape object
#' @return list
#' @seealso landscape.locus, landscape.locus.states
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   exampleland <- landscape.simulate(exampleland, 4)
#'   landscape.states(1,exampleland)
#'   rm(exampleland)
#' 
#' @export landscape.states
"landscape.states" <-
function(lnum=1,Rland)
  {
    if (is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          ain<-c();
          sta<-c();
          locin <- landscape.locus(lnum,Rland)[,c(-1:-landscape.democol())]
#          print(locin)
          ainds <- unique(c(locin))
#          print(ainds)
          for (i in 1:length(Rland$loci[[lnum]]$alleles))
            {
              if (Rland$loci[[lnum]]$alleles[[i]]$aindex %in% ainds)
                {
                  ain<-c(ain,Rland$loci[[lnum]]$alleles[[i]]$aindex);
                  sta<-c(sta,Rland$loci[[lnum]]$alleles[[i]]$state);
                }
            }
          sl <- list(aindex=ain,state=sta);
          if (landscape.ploidy(Rland)[lnum]>1)
            {
              st <- cbind(sl$state[match(locin[,1],sl$aindex)],sl$state[match(locin[,2],sl$aindex)])
            }
          else
            {
              st <- matrix(sl$state[match(c(locin),sl$aindex)],ncol=1)
            }
          cbind(landscape.locus(lnum,Rland)[,1:landscape.democol()],st)
        }
    
    
  }

