

library(kernelPop2)

rland <- NULL
  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland, h=50, s=2,np=0,totgen=20000)
  rland <- landscape.new.switchparam(rland,mp=0)
  rland <- landscape.new.floatparam(rland,s=0,seedscale=c(30,400),
                                    seedshape=c(1,400),seedmix=c(0.1),
                                    pollenscale=c(50,200),pollenshape=c(1,1),
                                    pollenmix=0.1 , asp=0.5)


  S <- matrix(c(
                  0.1, 0,
                0.25, 0
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

S <- matrix(0,ncol = (rland$intparam$habitats*rland$intparam$stages),
            nrow = (rland$intparam$habitats*rland$intparam$stages))

R <- S
M <- S
  
rights <- floor(seq(500,100000,length=50))
locs <- as.matrix(data.frame(lft=c(0,501,rights[c(-1:-2)]-diff(rights[-1])+1),
                   bot=rep(0,50),
                   rgt=rights,
                   top=10000))
  
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

  expmat <- cbind(c(0.5,0),c(0,0.5))
  hsq <- c(1,1)
  rland <- landscape.new.expression(rland,expmat=expmat,hsq=hsq)
rland <- landscape.new.gpmap(rland,
                             matrix(c(-1,0,1,
                                      -1,0,1,
                                      -1,0,1,
                                      1,-0.5,0.1,     #mixture
                                      -1,0,1),ncol=3,byrow=T),
                             matrix(c(-1,0,1,
                                      -1,0,1,
                                      0,-0.5,0.001    #reproduction
                                      ),ncol=3,byrow=T))
  initpopsize <- 500
rland <- landscape.new.individuals(rland,c(initpopsize,initpopsize,rep(0,98)))

###############################

l=rland
landscape.plot.locations(l)
for (i in 1:200)
{
    if (dim(l$individuals)[1]>0) l.old=l
    l=landscape.simulate(l,1)
    landscape.plot.locations(l)
    print(i)
    print(dim(l$individuals))
    print(landscape.allelefreq(l) )
    print(colMeans(landscape.phenotypes.c(l)))
}
