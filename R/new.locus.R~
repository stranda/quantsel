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
    as.integer(.Call("num_loci_poss",PACKAGE="kernelPop2"))
  }
