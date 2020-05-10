##get the density of individuals in each population (per sq grid unit)
##
#' @export
#'
popdens <- function(l)
{
    popsize <- c(table(landscape.populations(l)))
    poparea <- sapply(1:length(popsize),function(p)
    {
        (l$demography$epochs[[1]]$rightx[p]-l$demography$epochs[[1]]$leftx[p])*
            (l$demography$epochs[[1]]$topy[p]-l$demography$epochs[[1]]$boty[p])
        
    })
    popsize/poparea
}

