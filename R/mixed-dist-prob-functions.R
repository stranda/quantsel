#
# R functions to implement a mixed distribution with weibull and
# zero-truncated (positive only) gaussian
#

#density function
dmixed <- function(x,sc1=10,sh1=1,sc2=200,sh2=50,mix=0.75)
	{
			mix*dweibull(x,shape=sh1,scale=sc1)+(1-mix)*dnorm(x,mean=sc2,sd=sh2)/(1-pnorm(0,mean=sc2,sd=sh2))
	}

#cumulative probability function
#
pmixed <- function(x,sc1=10,sh1=1,sc2=200,sh2=50,mix=0.75)
	{	
			mix*pweibull(x,shape=sh1,scale=sc1)+(1-mix)*pnorm(x,mean=sc2,sd=sh2)/(1-pnorm(0,mean=sc2,sd=sh2))	}


#
# quantile function
# uses optimize.  Slow but works
qmixed <- function(p,sc1=10,sh1=1,sc2=200,sh2=50,mix=0.75)
	{
		f=function(q,scl1,shp1,scl2,shp2,mixp)
			{
				m=(100000*pmixed(q,sc1=scl1,sh1=shp1,sc2=scl2,sh2=shp2,mix=mixp)-100000*p)^2
				m
			}
		optimize(f,interval=c(0.0,max(3*sh2+sc2,5*sc1)),tol=1e-9,
			scl1=sc1,shp1=sh1,scl2=sc2,shp2=sh2,mixp=mix)$minimum
	}	

#
#Random number generating function
#
rmixed <- function(n=1,sc1=10,sh1=1,sc2=200,sh2=50,mix=0.75) {
  ifelse(runif(n)<=mix,
         rweibull(n,shape=sh1,scale=sc1),
         abs(rnorm(n,mean=sc2,sd=sh2))
         )
}




ldd <-function(scl=seq(100,2000,length=20),prop=0.1,mix=0.5)
{
	sapply(scl,function(x,prop,mix)
		{
			qmixed(p=0.99,sc2=x,sh2=prop*x,mix=mix)
		},prop=prop,mix=mix)
}


weibullvar <- function(scale=100,shape=2)
  {
    scale^2*(gamma(1+2/shape)-gamma(1+1/shape)^2)
  }
