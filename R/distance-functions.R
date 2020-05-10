#' @export
pollination.dist <- function(l)
  {
    sqrt((l$ind[,6]-l$ind[,8])^2+(l$ind[,7]-l$ind[,9])^2)
  }
#' @export
seed.dist <- function(l)
  {
    sqrt((l$ind[,6]-l$ind[,4])^2+(l$ind[,7]-l$ind[,5])^2)
  }
#' @export
midparent.dist <- function(l) #dispersal from the midparent location
{
    sqrt(
    (rowMeans(l$individuals[,c(6,8)])-l$individuals[,4])^2 +
    (rowMeans(l$individuals[,c(7,9)])-l$individuals[,5])^2
    )
}
#' @export
realized.selfing <- function(l)
  {
    pd <- pollination.dist(l)
    sum(pd==0)/length(pd)
  }
