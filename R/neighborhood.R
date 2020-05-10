##return the neighborhood size (both radius of 2*pi*vardisp and Nn)
##
#' @export
#'
landscape.neighborhood <- function(l)
{
    pd <- popdens(l)
    mid <- midparent.dist(l)
    seed <- seed.dist(l)
    pol <- pollination.dist(l)
    Nn <- 4 * pi * var(mid) * pd
    list(Nn=Nn,midmn=mean(mid),midsd=sd(mid),
         seedmn=mean(seed),seedsd=sd(seed),
         polmn=mean(pol),polsd=sd(pol))
}

