
"landscape.new.individuals" <-
function(rland, PopulationSizes)
  {
    if (!is.null(rland$intparam)&(rland$intparam$habitats*rland$intparam$stages==length(PopulationSizes)))
      rland <- .Call("populate_Rland",rland,PopulationSizes,PACKAGE="kernelPop2")
    else
      warning("either the Population Size vector does not match other landscape parameters or those parameters have not been declared")
    rland
  }

landscape.add.individuals <- function(rland, PopulationSizes)
  {
    if (!is.null(rland$intparam)&(rland$intparam$habitats*rland$intparam$stages==length(PopulationSizes)))
      {
        rland2 <- .Call("populate_Rland",rland,PopulationSizes,PACKAGE="kernelPop2")
        rland2$individuals[,3] <- rep(rland$intparam$currentgen,dim(rland2$individuals)[1])
      }
    else
      stop("either the Population Size vector does not match other landscape parameters or those parameters have not been declared")
    rland$individuals <- rbind(rland2$individuals,rland$individuals)
    rland$individuals <- rland$individuals[order(rland$individuals[,1]),]
    rland
  }

