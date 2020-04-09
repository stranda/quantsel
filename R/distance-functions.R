pollination.dist <- function(l)
  {
    sqrt((l$ind[,6]-l$ind[,8])^2+(l$ind[,7]-l$ind[,9])^2)
  }
seed.dist <- function(l)
  {
    sqrt((l$ind[,6]-l$ind[,4])^2+(l$ind[,7]-l$ind[,5])^2)
  }
realized.selfing <- function(l)
  {
    pd <- pollination.dist(l)
    sum(pd==0)/length(pd)
  }
