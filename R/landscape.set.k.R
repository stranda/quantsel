#' use a holosim spatiotemporal surface to set K in pops
#' @param N base level Ne that gets multiplied by per cell suit.
#' @param hsl landscape object used in holosimCell
#' @export
calc.k.vec <- function(N, hsl, gen, years=21000)
{    
    hs=hsl$hab_suit
    periods=dim(hs)[1]
    tbl=as.numeric((cut(1:years,breaks=periods)))
    k=N * hs[tbl[gen],]
    k[is.na(k)]=0
    k
} #rnd function
