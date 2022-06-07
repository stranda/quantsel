#' calculate summary statistics on habitats in 2d landscape
#' @param l quantsel landscape
#' @param e epoch (default 1, usually fine)
#' @return a named list of cell stats
#' @export
landscape.hab.stats <- function(l, e=1)
{
    if (is.landscape(l))
    {
       leftx=l$demography$epochs[[e]]$leftx
       rightx=l$demography$epochs[[e]]$rightx
       topy=l$demography$epochs[[e]]$topy
       boty=l$demography$epochs[[e]]$boty
       list(numhabs=length(leftx),
            mnareas = mean((rightx-leftx)*(topy-boty)),
            mnx = mean(rightx-leftx),
            mny = mean(topy-boty),
            mndiag=mean(sqrt((rightx-leftx)^2+(topy-boty)^2))
            )
    } else {NULL}
}

