#
#mixed weibull and normal distribution densities
#
#
dmixed <- function(x,scale,mix,scale2,shape2,log=F,order=F)
  {
    res <-  dweibull(abs(x),scale=scale,shape=1)*mix + dnorm(abs(x),mean=scale2,sd=shape2)*(1-mix)
    if (order)
      if (scale>scale2)
        res <- 0
    if (log)
      {
        res <- log(res)
      }
    res
  }
#
#mixed weibull and normal distribution cumulative densities
#
#
pmixed <- function(x,scale,mix,scale2,shape2,log=F,order=F)
  {
    res <-  pweibull(abs(x),scale=scale,shape=1)*mix + pnorm(abs(x),mean=scale2,sd=shape2)*(1-mix)
    if (order)
      if (scale>scale2)
        res <- 0
    if (log)
      {
        res <- log(res)
      }
    res
  }

#
#mixed weibull and normal distribution quantile
#
#
qmixed.sc2 <- function(p,scale,mix,scale2)
  {
    sapply(scale2,function(scl2,p,scale,mix)
           {
             as.numeric(optimize(f=function(x,p=p,scale,mix,scale2,shape2)
                                 {
                                   (pmixed(x,scale,mix,scale2,shape2)-p)^2
                                 },interval=c(0,10000),
                                 p=p,scale=scale,mix=mix,scale2=scl2,
                                 shape2=0.8*scl2)$minimum)
           },
           p=p,scale=scale,mix=mix)
  }
