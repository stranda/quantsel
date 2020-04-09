"landscape.new.example" <- function()
{
  rland <- NULL
  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland, h=20, s=2,np=0,totgen=20000)
  rland <- landscape.new.switchparam(rland,mp=0)
  rland <- landscape.new.floatparam(rland,s=0,seedscale=c(90,600),
                                    seedshape=c(1,200),seedmix=c(1),
                                    pollenscale=c(50,50),pollenshape=c(1,1),
                                    pollenmix=1 )


  S <- matrix(c(
                  0.33, 0,
                0.0175, 0
                ), byrow=T, nrow = 2)
  R <- matrix(c(
                0, 400,
                0,   0
                ), byrow=T, nrow = 2)
  M <- matrix(c(
                0, 0,
                0, 1
                ), byrow=T, nrow = 2)

  rland <- landscape.new.local.demo(rland,S,R,M)
  
  S <- matrix(rep(0,(2*rland$intparam$habitat)^2), nrow = 2*rland$intparam$habitat)
  R <- matrix(rep(0,(2*rland$intparam$habitat)^2), nrow = 2*rland$intparam$habitat)
  M <- matrix(rep(0,(2*rland$intparam$habitat)^2), nrow = 2*rland$intparam$habitat)

  locs <- landscape.generate.locations(npop=rland$intparam$habitat,
                             xrange=c(0,15000),yrange=c(0,15000),
                             sizexkernel=c(500,65),sizeykernel=c(500,65)
                             )
  
  rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                     carry=(1.50 * (sqrt((locs[,3]-locs[,1])*(locs[,4]-locs[,2])))),
                     extinct=rep(0.05,rland$intparam$habitat),
                     leftx=locs[,1],
                     rightx=locs[,3],
                     boty=locs[,2],
                     topy=locs[,4],
                     maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))
  

  rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)
  rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)

  expmat <- cbind(c(1,0),c(0,1))
  hsq <- c(1,1)
  rland <- landscape.new.expression(rland,expmat=expmat,hsq=hsq)
  rland <- landscape.new.gpmap(rland)
  initpopsize <- 150
  rland <- landscape.new.individuals(rland,
                         round(runif(2*rland$intparam$habitat,
                                     min=0,
                                     max=initpopsize/1.5)))
  rland
}
