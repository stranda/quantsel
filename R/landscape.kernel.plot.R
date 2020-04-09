#
# plots a contour diagram of the 2d kernel given a landscape/currently non-functional
#
landscape.kernel.plot <- function(rland)
  {
    xlim <- (max(rland$demography$epochs[[1]]$leftx)-min(rland$demography$epochs[[1]]$rightx))/2
    ylim <- (max(rland$demography$epochs[[1]]$topy)-min(rland$demography$epochs[[1]]$boty))/2
    zmat <- outer(seq(-xlim,xlim,length=100),
                  seq(-xlim,xlim,length=100),
                  function(x,y,mu,mix,mu2,shape2)
                  {
                    log(dmixed(abs(x),mu,mix,mu2,shape2))
                  },
                  mu=rland$floatparam$seedmu,
                  mix=rland$floatparam$seedmix,
                  mu2=rland$floatparam$seedmu2,
                  shape2=rland$floatparam$seedshape2
                )
    contour(xlim,ylim,zmat)
  }
