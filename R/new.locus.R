"landscape.new.locus.old" <-
function(rland,type=0,ploidy=1,mutationrate=0,transmission=1,numalleles=2,frequencies=NULL,states=NULL)
{
  if (!(is.list(rland$loci)))
    {
      rland$loci <- list(list(type=0,ploidy=0,trans=0,rate=0,alleles=0))
      locusnum <- 1
    }
  else
    {
      locusnum <- length(rland$loci) + 1
      rland$loci[[locusnum]] <- list(type=0,ploidy=0,rate=0,trans=0,alleles=0)
    }

  rland$intparam$locusnum <- locusnum
  
  if(type >= 0 && type <= 1)
    {
      rland$loci[[locusnum]]$type <- typelookup(type)
    }
  else
    {
      stop("Invalid type of locus")
    }

  if(ploidy == 1 || ploidy == 2)
    {
      rland$loci[[locusnum]]$ploidy <- ploidy
    }
  else
    {
      stop("Invalid ploidy count")
    }

  rland$loci[[locusnum]]$rate <- mutationrate

  if(transmission == 0 || transmission == 1)
    {
      rland$loci[[locusnum]]$trans <- as.integer(transmission)
    }
  else
    {
      stop("Invalid transmission number")
    }

  if(numalleles >= 0)
    {
      rland$loci[[locusnum]]$alleles <- makealleles(type,numalleles,allelesize,frequencies,states)
    }
  else
    {
      stop("Need non-negative numbers of alleles")
    }

  rland
}





#' Add a locus
#' 
#' Add a locus to a Rmetasim landscape object
#' 
#' 
#' @param rland partially created landscape object, required
#' @param type (default=0) type of locus, 0=Infinite Allele mutation model
#' (Integer), 1=Stepwise mutation model (Integer) state, 2=DNA base (variable
#' length string state)
#' @param ploidy (default=1) locus ploidy, 1 or 2
#' @param mutationrate (default=0) probability of mutation per generation, less
#' than or equal to 1
#' @param transmission (default=1) 1=uniparental inheritance, 0=biparental
#' inheritance
#' @param numalleles (default=2) number of different alleles at the time of
#' creation
#' @param allelesize (default=50) length of DNA strings if type=2
#' @param frequencies (default=NULL) vector of frequencies for each allele,
#' must be numalleles long and add up to 1, if NULL frequencies are equally
#' distributed
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.empty()
#'   exampleland <- landscape.new.intparam(exampleland, s=2, h=2)
#'   exampleland <- landscape.new.floatparam(exampleland)
#'   exampleland <- landscape.new.switchparam(exampleland)
#' 
#'   exampleland <- landscape.new.locus(exampleland,type=2,ploidy=2,mutationrate=.001,numalleles=5,allelesize=100)
#' 
#'   exampleland$loci
#' 
#'   rm(exampleland)
#' 
#' @export landscape.new.locus
landscape.new.locus <- function (rland, type = 0, ploidy = 1, mutationrate = 0, transmission = 1, 
    numalleles = 2, allelesize = 50, frequencies = NULL, states = NULL) 
{
    if (!(is.list(rland$loci))) {
        rland$loci <- list(list(type = 0, ploidy = 0, trans = 0, 
            rate = 0, alleles = 0))
        locusnum <- 1
    }
    else {
        locusnum <- length(rland$loci) + 1
        if (landscape.maxloci()<locusnum)
           stop(paste0("Trying to create more loci(",locusnum,")than currently allowed by kernelPop:(",landscape.maxloci(),")\nTo fix either: reduce number of loci\nor change the constant MAXLOCI in 'src/const.h' to a value slightly greater than\n\tyou need and recompile rmetasim.  If recompiling is a problem, please contact author"))
        rland$loci[[locusnum]] <- list(type = 0, ploidy = 0, 
            rate = 0, trans = 0, alleles = 0)
    }
    rland$intparam$locusnum <- locusnum
    if (type >= 0 && type <= 2) {
        rland$loci[[locusnum]]$type <- as.integer(typelookup(type))
    }
    else {
        stop("Invalid type of locus")
    }
    if (ploidy == 1 || ploidy == 2) {
        rland$loci[[locusnum]]$ploidy <- as.integer(ploidy)
    }
    else {
        stop("Invalid ploidy count")
    }
    rland$loci[[locusnum]]$rate <- mutationrate
    if (transmission == 0 || transmission == 1) {
        rland$loci[[locusnum]]$trans <- as.integer(transmission)
    }
    else {
        stop("Invalid transmission number")
    }
    if (numalleles >= 0) {
        rland$loci[[locusnum]]$alleles <- makealleles(type, numalleles, 
            allelesize, frequencies, states)
    }
    else {
        stop("Need non-negative numbers of alleles")
    }
    rland
}



#
# a convenience function that provides the highest number of loci possible (value of MAXLOCI in 'const.h')column number for demographic information in
# the individuals matrix
#
landscape.maxloci <- function()
  {
    as.integer(.Call("num_loci_poss",PACKAGE="quantsel"))
  }
