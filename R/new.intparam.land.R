#' Create a set of integer parameters
#' 
#' Create a set of integer parameters for a Rmetasim landscape.
#' 
#' 
#' @param rland skeletion of landscape object, required
#' @param h habitats (default=1), the number of different subpopulations within
#' the landscape
#' @param s stages (default=1), the number of stages in the life cycle of the
#' organism
#' @param cg currentgen (default=0), the current generation the simulation has
#' reached
#' @param ce currentepoch (default=0), the current epoch the simulation has
#' reached
#' @param totgen totoalgens (default=1000), the total number of generations to
#' simulate
#' @param maxland maxlandsize(default=200000), the maxium number of individuals
#' that can exist in the simulation
#' @keywords misc
#' @examples
#' 
#'   ## Defaults
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.intparam(exampleland)
#'   exampleland$intparam
#' 
#'   ## 2 habitats, 3 stage lifecycle, 1000000 generations, maximum 1000000 individuals
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.intparam(exampleland,h=2,s=2,totgen=1000000,maxland=1000000)
#'   exampleland$intparam
#' 
#'   rm(exampleland)
#' 
#' @export landscape.new.intparam
"landscape.new.intparam" <-
function(rland,h=1,s=1,cg=0,ce=0,totgen=1000,maxland=300000,np=0,rows=0,cols=0)
{
  rland$intparam <- list(h,s,0,0,cg,ce,totgen,0,maxland,np,rows,cols)
  names(rland$intparam) <- c("habitats","stages","locusnum","numepochs","currentgen","currentepoch","totalgens","numdemos","maxlandsize","nphen","rows","cols")
  rland
}

