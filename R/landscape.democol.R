#
# a convenience function that provides the highest column number for demographic information in
# the individuals matrix
#




#' return largest demographic column from a landscape
#' 
#' return largest demographic column from a landscape
#' 
#' Useful to write functions that will be insensitive to some changes in the
#' individuals object (mainly addition of non-genetic information)
#' 
#' @return a scalar integer representing the largest column of demographic
#' information in a landscape's individuals object
#' @seealso landscape.locus
#' @keywords misc
#' @export landscape.democol
landscape.democol <- function()
{
    as.integer(9)
  }
