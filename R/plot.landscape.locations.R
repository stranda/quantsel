#
#
#adults selects plotting of adults only
#




#' plots the locations of all habitats and individuals
#' 
#' 
#' Plots the landscape.  Gives different habitats different colors, though the
#' colors cycle through rapidly when the number of habitats is large.
#' Offspring always have the color of their \emph{mother's} habitat.  This
#' provides a quick way to assess inter-habitat movement visually
#' 
#' 
#' @param rland landscape
#' @param adults vector of stages to select and plot, if NULL, plot all stages
#' @param label boolean to indicate if cells labeled (default FALSE)
#' @return NULL
#' @keywords misc
#' @export landscape.plot.locations
landscape.plot.locations <- function(rland,adults=c(NULL),label=FALSE)
  {
    if (is.landscape(rland,FALSE))
      {
        plot(1,1,type="n",
             xlim=c(min(rland$demography$epochs[[1]]$leftx),max(rland$demography$epochs[[1]]$rightx)),
             ylim=c(min(rland$demography$epochs[[1]]$boty),max(rland$demography$epochs[[1]]$topy)),
             xlab="X coordinate",ylab="Y coordinate",
             main=paste("landscape state at generation",rland$intparam$currentgen))
        for (i in 1:rland$intparam$habitats)
          {
            rect(rland$demography$epochs[[1]]$leftx[i],
                 rland$demography$epochs[[1]]$boty[i],
                 rland$demography$epochs[[1]]$rightx[i],
                 rland$demography$epochs[[1]]$topy[i],
                 lwd=2,border=i)
            if (label) text(x=(rland$demography$epochs[[1]]$rightx[i]+rland$demography$epochs[[1]]$leftx[i])/2,
                 y=(rland$demography$epochs[[1]]$boty[i]+rland$demography$epochs[[1]]$topy[i])/2,
                 as.character(i-1),cex=0.5)
          }
        if (length(landscape.populations(rland))>1)
          {
            if (!is.null(adults[1]))
              {
                rland$individuals <- rland$individuals[rland$individuals[,1] %in% adults,]
              }
            icol <- getpopulations(rland,rland$individuals[,c(6)],rland$individuals[,c(7)])
            points(rland$individuals[,c(4,5)],type="p",
                   pch=0+(rland$individuals[,1] - (rland$intparam$stages*(landscape.populations(rland)-1))),
                   col=icol+1,
                   cex=0.5)
          }
      }
    else
      {
        print ("no landscape to plot")
      }
}

                                        #Take a vector of x,y and return the population numbers
getpopulations <- function(rland,x,y)
{
  popmat <- sapply(1:rland$intparam$habitats,function(i,x,y)
               {
                 as.logical((rland$demography$epochs[[1]]$leftx[i]<=x)*(rland$demography$epochs[[1]]$rightx[i]>=x)*(rland$demography$epochs[[1]]$boty[i]<y)*(rland$demography$epochs[[1]]$topy[i]>y))
               },x=x,y=y)

  rv <- unlist(apply(popmat,1,which))
  rv[(x==0)&(y==0)] <- NA
  rv
}
                                        #Take a vector of x,y and return the population numbers
getpopulations.old <- function(rland,x,y)
{
  
  popmat <- sapply(1:rland$intparam$habitats,function(i,x,y)
               {
                 as.logical((rland$demography$epochs[[1]]$leftx[i]<=x)*(rland$demography$epochs[[1]]$rightx[i]>=x)*(rland$demography$epochs[[1]]$boty[i]<y)*(rland$demography$epochs[[1]]$topy[i]>y))
               },x=x,y=y)
  retval <- unlist(apply(popmat,1,which))[unlist(apply(popmat,1,which))>-1]
  retval[sapply((unlist(apply(popmat,1,which))[unlist(apply(popmat,1,which))>-1]),is.null)] <- 0
  unlist(retval)

}
